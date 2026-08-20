#!/usr/bin/env python3
"""RNA-seq quantification with pytximport.

Replaces deprecated_Quantification.qmd (kept in the repo as reference).
Merged into the bash pipeline as "Phase 4: Quantification" of
run_rnaseq.sh, but can also be run standalone from the project root:

    python quantification.py config.yml

It imports Salmon quant.sf files with pytximport, transfers transcript-level
counts/TPM to gene-level, optionally adds featureCounts (BAM) counts/TPM, and
writes date-stamped xlsx files into the quantification directory:

    <date>_transcript_expression.xlsx      (sheets: tpm, counts)
    <date>_gene_expression.xlsx            (sheets: tpm, counts)
    <date>_gene_expression_bam.xlsx        (only if featureCounts output exists)

The sheet layout matches the old R output (columns: ensembl_id, symbol, type,
then one column per sample in samplesheet order), so
DESeq2.qmd and Enrichment.qmd read these files unchanged.

Idempotent: skipped when today's <date>_gene_expression.xlsx already exists.
"""

import argparse
import gzip
import logging
import os
import re
import sys
from datetime import date

import numpy as np
import pandas as pd

from pytximport import tximport

log = logging.getLogger("quantification")

# Species -> (index subdir, gtf file name). Mirrors helper.sh setup_genome().
SPECIES_GTF = {
    "chm13": ("hs/chm13", "chm13v2.0.gtf.gz"),
    "hs": ("hs/v49", "gencode.v49.basic.annotation.gtf.gz"),
    "mm": ("mm/vM38", "gencode.vM38.basic.annotation.gtf.gz"),
    "fly": ("fly", "Drosophila_melanogaster.BDGP6.54.63.gtf.gz"),
    "cel": ("cel", "Caenorhabditis_elegans.WBcel235.63.gtf.gz"),
}

VERSION_RE = re.compile(r"\.[0-9]+$")


# ============================= Config =====================================
def resolve_path(root_dir, path):
    """Join a (possibly absolute) config path onto root_dir like the bash pipeline."""
    if os.path.isabs(path):
        return path
    return os.path.join(root_dir, path)


def parse_config(config_path):
    import yaml

    with open(config_path) as fh:
        cfg = yaml.safe_load(fh) or {}

    root_dir = cfg.get("root_dir", ".")
    outdir = cfg.get("outdir", "result")
    species = cfg.get("species", "hs")

    paths = {
        "root_dir": root_dir,
        "samplesheet": resolve_path(root_dir, cfg.get("samplesheet", "samplesheet.csv")),
        "salmon_dir": resolve_path(root_dir, cfg.get("salmon_dir", os.path.join(outdir, "03_salmon"))),
        "bam_count_dir": resolve_path(root_dir, cfg.get("bam_count_dir", os.path.join(outdir, "05_featurecounts"))),
        "quantification_dir": resolve_path(root_dir, cfg.get("quantification_dir", "quantification")),
    }

    if "gtf_file" in cfg:
        paths["gtf_file"] = resolve_path(root_dir, cfg["gtf_file"])
    else:
        index_rootdir = cfg.get("index_rootdir", "/mnt/d/index")
        if species not in SPECIES_GTF:
            log.error("Unknown species '%s' (set gtf_file in config.yml)", species)
            sys.exit(1)
        subdir, gtf_name = SPECIES_GTF[species]
        paths["gtf_file"] = os.path.join(index_rootdir, subdir, gtf_name)

    return paths


# ============================= GTF parsing ================================
_ATTR_VALUE_RE = {
    key: re.compile(r'%s\s+"([^"]*)"' % key)
    for key in (
        "gene_id",
        "gene_name",
        "gene_type",
        "gene_biotype",
        "transcript_id",
        "transcript_name",
        "transcript_type",
        "transcript_biotype",
    )
}


def _attr(attr_str, key):
    m = _ATTR_VALUE_RE[key].search(attr_str)
    return m.group(1) if m else ""


def strip_version(gene_id):
    return VERSION_RE.sub("", gene_id)


def parse_gtf(gtf_file):
    """Return (transcript_gene_map, gene_annotation, transcript_annotation).

    gene_annotation / transcript_annotation have columns:
        ensembl_id, symbol, type
    and are distinct on ensembl_id (mirrors the R 'pharse_gtf' output).
    """
    opener = gzip.open if gtf_file.endswith(".gz") else open

    gene_rows = {}
    tx_rows = {}
    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            attr = fields[8]
            if feature == "gene":
                gid = strip_version(_attr(attr, "gene_id"))
                if not gid:
                    continue
                gene_rows.setdefault(gid, {
                    "ensembl_id": gid,
                    "symbol": _attr(attr, "gene_name"),
                    "type": _attr(attr, "gene_type") or _attr(attr, "gene_biotype"),
                })
            elif feature == "transcript":
                tid = strip_version(_attr(attr, "transcript_id"))
                gid = strip_version(_attr(attr, "gene_id"))
                if not tid or not gid:
                    continue
                tx_rows.setdefault(tid, {
                    "ensembl_id": tid,
                    "symbol": _attr(attr, "transcript_name"),
                    "type": _attr(attr, "transcript_type") or _attr(attr, "transcript_biotype"),
                    "gene_id": gid,
                })

    t2g = pd.DataFrame(
        [{"transcript_id": t, "gene_id": r["gene_id"]} for t, r in tx_rows.items()]
    ).drop_duplicates("transcript_id")

    gene_annotation = pd.DataFrame(list(gene_rows.values()))
    transcript_annotation = pd.DataFrame(
        [{k: v for k, v in r.items() if k != "gene_id"} for r in tx_rows.values()]
    )
    return t2g, gene_annotation, transcript_annotation


# ============================= Samplesheet ================================
def parse_samplesheet(path):
    sample = pd.read_csv(path)
    if "sample" not in sample.columns:
        log.error("samplesheet '%s' missing 'sample' column", path)
        sys.exit(1)
    samples = [str(x) for x in sample["sample"].tolist()]
    if len(set(samples)) != len(samples):
        log.error("duplicated sample names in samplesheet")
        sys.exit(1)
    return samples


# ============================= Quantification =============================
def matrix_from_anndata(ad, key, samples):
    """Extract a features x samples DataFrame from pytximport AnnData output."""
    if key == "counts":
        if "counts" in ad.obsm:
            m = np.asarray(ad.obsm["counts"])
        elif ad.X is not None:
            m = np.asarray(ad.X)
        else:
            raise KeyError(f"'counts' not in pytximport output (obsm keys: {sorted(ad.obsm)})")
    elif key in ad.obsm:
        m = np.asarray(ad.obsm[key])
    else:
        raise KeyError(f"key '{key}' not in pytximport output (obsm keys: {sorted(ad.obsm)})")
    if m.ndim == 1:
        m = m.reshape(1, -1)
    if m.shape[0] != len(samples):
        raise ValueError("sample count mismatch in pytximport output")
    return pd.DataFrame(m.T, index=list(ad.var_names), columns=samples)


def annotate_expression(expr, annotation, decimals):
    """Join expression matrix to annotation, keep expr row order (R right_join semantics).

    expr: DataFrame with index = ensembl ids (possibly versioned), columns = samples.
    Returns DataFrame with columns: ensembl_id, symbol, type, <samples>.
    """
    df = expr.rename_axis("ensembl_id").reset_index()
    df["ensembl_id"] = df["ensembl_id"].astype(str).map(strip_version)
    sample_cols = [c for c in df.columns if c != "ensembl_id"]
    out = df.merge(annotation, how="left", on="ensembl_id")
    out = out[["ensembl_id", "symbol", "type"] + sample_cols]
    out[sample_cols] = out[sample_cols].round(decimals)
    return out


def write_xlsx(path, tpm, counts):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        tpm.to_excel(writer, sheet_name="tpm", index=False)
        counts.to_excel(writer, sheet_name="counts", index=False)


def run_salmon_quantification(paths, samples, t2g, annotations, today):
    quant_files = [os.path.join(paths["salmon_dir"], s, "quant.sf") for s in samples]
    missing = [f for f in quant_files if not os.path.isfile(f)]
    if missing:
        log.error("sample(s) in samplesheet are not identical with salmon output:")
        for f in missing:
            log.error("  missing: %s", f)
        sys.exit(1)

    log.info("Importing %d Salmon quantification files", len(quant_files))
    tx_txi = tximport(
        quant_files, data_type="salmon", transcript_gene_map=t2g,
        output_type="anndata", counts_from_abundance=None, return_transcript_data=True,
    )
    gene_txi = tximport(
        quant_files, data_type="salmon", transcript_gene_map=t2g,
        output_type="anndata", counts_from_abundance=None,
    )

    gene_annotation, transcript_annotation = annotations

    log.info("Writing transcript-level expression")
    tx_counts = matrix_from_anndata(tx_txi, "counts", samples)
    tx_tpm = matrix_from_anndata(tx_txi, "abundance", samples)
    write_xlsx(
        os.path.join(paths["quantification_dir"], f"{today}_transcript_expression.xlsx"),
        annotate_expression(tx_tpm, transcript_annotation, 2),
        annotate_expression(tx_counts, transcript_annotation, 0),
    )

    log.info("Summarizing to gene-level expression")
    gene_counts = matrix_from_anndata(gene_txi, "counts", samples)
    gene_tpm = matrix_from_anndata(gene_txi, "abundance", samples)
    write_xlsx(
        os.path.join(paths["quantification_dir"], f"{today}_gene_expression.xlsx"),
        annotate_expression(gene_tpm, gene_annotation, 2),
        annotate_expression(gene_counts, gene_annotation, 0),
    )

    return True


def run_featurecounts_quantification(paths, samples, gene_annotation, today):
    counts_file = os.path.join(paths["bam_count_dir"], "gene_counts.txt")
    if not os.path.isfile(counts_file):
        return False

    log.info("Reading featureCounts output: %s", counts_file)
    fc = pd.read_csv(counts_file, sep="\t", comment="#")
    fc = fc.drop(columns=[c for c in ("Chr", "Start", "End", "Strand") if c in fc.columns])

    rename = {}
    for c in fc.columns:
        if c not in ("Geneid", "Length") and c.endswith(".bam"):
            rename[c] = os.path.basename(c)[: -len(".bam")]
    fc = fc.rename(columns=rename)

    cols = [c for c in ("Geneid", "Length") if c in fc.columns] + [
        s for s in samples if s in fc.columns
    ]
    fc = fc[cols].rename(columns={"Geneid": "ensembl_id"})
    fc["ensembl_id"] = fc["ensembl_id"].astype(str).map(strip_version)

    sample_cols = [s for s in samples if s in fc.columns]
    length = fc["Length"].astype(float)

    # counts + annotation (R right_join semantics)
    counts = fc[["ensembl_id"] + sample_cols].copy()
    counts = counts.merge(gene_annotation, how="left", on="ensembl_id")
    counts = counts[["ensembl_id", "symbol", "type"] + sample_cols]
    counts[sample_cols] = counts[sample_cols].round(0)

    # TPM from counts / gene length
    rates = fc[sample_cols].div(length, axis=0)
    tpm = rates.div(rates.sum(axis=0), axis=1).mul(1e6)
    tpm = tpm.assign(ensembl_id=fc["ensembl_id"].values)
    tpm = tpm.merge(gene_annotation, how="left", on="ensembl_id")
    tpm = tpm[["ensembl_id", "symbol", "type"] + sample_cols]
    tpm[sample_cols] = tpm[sample_cols].round(2)

    write_xlsx(
        os.path.join(paths["quantification_dir"], f"{today}_gene_expression_bam.xlsx"),
        tpm, counts,
    )
    return True


# ============================= Main =======================================
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", default="config.yml",
                        help="config.yml path (default: config.yml)")
    parser.add_argument("--date", default=None, help="override output date prefix (YYYY-MM-DD)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, stream=sys.stderr,
        format="[%(asctime)s] [INFO] %(message)s", datefmt="%Y-%m-%d %H:%M",
    )

    if not os.path.isfile(args.config):
        log.error("config file not found: %s", args.config)
        sys.exit(1)

    paths = parse_config(args.config)
    today = args.date or date.today().isoformat()

    if not os.path.isfile(paths["samplesheet"]):
        log.error("samplesheet not found: %s", paths["samplesheet"])
        sys.exit(1)
    samples = parse_samplesheet(paths["samplesheet"])

    gene_out = os.path.join(paths["quantification_dir"], f"{today}_gene_expression.xlsx")
    if os.path.isfile(gene_out):
        log.info("Quantification output already exists, skipping: %s", gene_out)
        return

    log.info("Parsing annotation: %s", paths["gtf_file"])
    if not os.path.isfile(paths["gtf_file"]):
        log.error("GTF file not found: %s", paths["gtf_file"])
        sys.exit(1)
    t2g, gene_annotation, transcript_annotation = parse_gtf(paths["gtf_file"])
    log.info("  genes: %d, transcripts: %d", len(gene_annotation), len(transcript_annotation))

    ran_salmon = False
    if os.path.isdir(paths["salmon_dir"]):
        ran_salmon = run_salmon_quantification(
            paths, samples, t2g, (gene_annotation, transcript_annotation), today)
    else:
        log.info("no salmon output directory found, skip it: %s", paths["salmon_dir"])

    ran_bam = False
    if os.path.isdir(paths["bam_count_dir"]):
        ran_bam = run_featurecounts_quantification(paths, samples, gene_annotation, today)
    else:
        log.info("no featurecounts directory found, skip it: %s", paths["bam_count_dir"])

    if not (ran_salmon or ran_bam):
        log.error("nothing to quantify (no salmon or featureCounts output)")
        sys.exit(1)

    log.info("Quantification completed -> %s", paths["quantification_dir"])


if __name__ == "__main__":
    main()