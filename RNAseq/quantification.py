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

Idempotent per product: each output xlsx is written only when missing, so a
partial day (e.g. Salmon done, featureCounts added later) can be completed by
re-running without redoing the earlier work.
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

# pytximport counts_from_abundance mode for the GENE-level matrix:
#   None                -> per-gene counts are simple sums of transcript counts.
#                          DESeq2 then REQUIRES length-based offsets (e.g. via
#                          DESeqDataSetFromTximport); plain summed counts in a
#                          bare DESeqDataSetFromMatrix bias against genes with
#                          preferentially long transcripts.
#   "length_scaled_tpm" -> tximport-recommended mode when the count matrix is
#                          fed to DESeq2 as-is (no offsets): counts are scaled
#                          by average transcript length so a plain matrix is
#                          method-correct.
# The xlsx consumers here read the counts sheet as a plain matrix, so
# length-scaled TPM is the method-correct default for the gene-level sheet.
# The transcript-level sheet keeps raw summed counts.
COUNTS_FROM_ABUNDANCE = "length_scaled_tpm"

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
    # dtype=str: pandas type inference turns a sample named "01" into 1,
    # which then fails to match the salmon output directory
    sample = pd.read_csv(path, dtype=str)
    if "sample" not in sample.columns:
        log.error("samplesheet '%s' missing 'sample' column", path)
        sys.exit(1)
    samples = [str(x) for x in sample["sample"].dropna().tolist()]
    if len(set(samples)) != len(samples):
        log.error("duplicated sample names in samplesheet")
        sys.exit(1)
    return samples


# ============================= Quantification =============================
def check_sample_order(ad, samples):
    """Fail loudly when pytximport's row order does not match the samplesheet.

    The expression matrix is built by zipping pytximport rows with the
    samplesheet order; if the two ever disagree (e.g. a pytximport version
    starts sorting its output), samples would be silently mislabelled —
    the most dangerous failure mode here.
    """
    obs = [str(x) for x in ad.obs_names]
    base = [os.path.basename(o).removesuffix(".sf") for o in obs]
    if set(base) == set(samples):
        if base != samples:
            raise ValueError(
                f"pytximport sample order {base} does not match samplesheet order {samples}"
            )
    else:
        # obs names are not resolvable to sample names (e.g. all "quant.sf"):
        # fall back to the count-only check in matrix_from_anndata
        log.debug("cannot map pytximport obs names to samples; skipping order check: %s", obs[:3])
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
    check_sample_order(tx_txi, samples)
    gene_txi = tximport(
        quant_files, data_type="salmon", transcript_gene_map=t2g,
        output_type="anndata", counts_from_abundance=COUNTS_FROM_ABUNDANCE,
    )
    check_sample_order(gene_txi, samples)

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
    missing = [s for s in samples if s not in fc.columns]
    if missing:
        log.warning(
            "featureCounts output is missing sample(s) %s; they will be "
            "absent from the bam xlsx", missing,
        )
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
    bam_out = os.path.join(paths["quantification_dir"], f"{today}_gene_expression_bam.xlsx")

    # Idempotency is PER PRODUCT: an early return keyed on the salmon output
    # alone used to make the featureCounts xlsx unreachable when Salmon had
    # already run the same day.
    skip_salmon = os.path.isfile(gene_out)
    skip_bam = os.path.isfile(bam_out)
    if skip_salmon:
        log.info("Gene-level Salmon output already exists, skipping salmon import: %s", gene_out)
    if skip_bam:
        log.info("featureCounts output already exists, skipping: %s", bam_out)
    if skip_salmon and skip_bam:
        log.info("All quantification outputs already exist for %s", today)
        return

    log.info("Parsing annotation: %s", paths["gtf_file"])
    if not os.path.isfile(paths["gtf_file"]):
        log.error("GTF file not found: %s", paths["gtf_file"])
        sys.exit(1)
    t2g, gene_annotation, transcript_annotation = parse_gtf(paths["gtf_file"])
    log.info("  genes: %d, transcripts: %d", len(gene_annotation), len(transcript_annotation))

    ran_salmon = skip_salmon  # existing output counts as "done"
    if not skip_salmon:
        if os.path.isdir(paths["salmon_dir"]):
            ran_salmon = run_salmon_quantification(
                paths, samples, t2g, (gene_annotation, transcript_annotation), today)
        else:
            log.info("no salmon output directory found, skip it: %s", paths["salmon_dir"])

    ran_bam = skip_bam
    if not skip_bam:
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