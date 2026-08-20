#!/usr/bin/env python3
"""plot_pca.py - PCA of multiBigwigSummary/multiBamSummary output.

Replacement for deepTools plotPCA (which in 3.5.6 ignores --log2 and draws
zero-lines that squeeze the data cloud). This script:

  - reads the .npz from multiBigwigSummary/multiBamSummary (bins x samples)
  - drops zero-variance rows, keeps the top --ntop rows by variance
  - z-scores each sample column (mean 0, unit variance)
  - SVD -> eigenvalues = s^2, % variance explained
  - plots PC1 vs PC2 with adaptive axis limits, colored by group
    (from the samplesheet 'group' column; falls back to the sample-name
    prefix for samples missing from the sheet)
  - writes a .tab with per-sample loadings and eigenvalues

Usage:
  python plot_pca.py -in corData.npz -o pca.pdf --outFileNameData pca.tab [--ntop N] [--plotTitle T] [--samplesheet samplesheet.csv]
"""

import argparse
import csv
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt


def parse_arguments():
    p = argparse.ArgumentParser(
        description="PCA plot from multiBigwigSummary/multiBamSummary output",
        usage="python plot_pca.py -in coverage.npz -o pca.pdf [--ntop N] [--outFileNameData file.tab] [--samplesheet samplesheet.csv]",
    )
    p.add_argument("-in", "--corData", required=True, help="Coverage file (.npz) from multiBigwigSummary/multiBamSummary")
    p.add_argument("-o", "--plotFile", required=True, help="Output plot file (pdf/png/svg/eps)")
    p.add_argument("--outFileNameData", metavar="FILE", help="Write per-sample loadings and eigenvalues (Component, samples..., Eigenvalue)")
    p.add_argument("--ntop", type=int, default=100000, help="Use only the top N most variable rows. 0 = all rows. (Default: 100000)")
    p.add_argument("--samplesheet", metavar="FILE", help="CSV with sample and group columns (used for coloring)")
    p.add_argument("--plotTitle", default="", help="Title of the plot")
    p.add_argument("--plotHeight", type=float, default=10, help="Plot height in cm")
    p.add_argument("--plotWidth", type=float, default=10, help="Plot width in cm")
    return p.parse_args()


def group_prefix(label):
    m = re.match(r"^[^0-9-]+", str(label))
    return m.group(0) if m else str(label)


def load_sample_groups(samplesheet):
    """Return {sample: group} from the samplesheet CSV (handles CRLF)."""
    mapping = {}
    with open(samplesheet, newline="") as f:
        for row in csv.DictReader(f):
            sample = (row.get("sample") or "").strip()
            group = (row.get("group") or "").strip()
            if sample:
                mapping[sample] = group
    return mapping


def categorical_colors(n):
    """Categorical colors that auto-expand when there are many groups."""
    tab10 = plt.cm.tab10.colors[:10]
    if n <= 10:
        return list(tab10[:n])
    tab20 = plt.cm.tab20.colors[:20]
    if n <= 20:
        return list(tab20[:n])
    # >20 groups: golden-angle HSV sampling for maximally distinct hues
    phi = (1 + 5**0.5) / 2
    hues = [(i * phi) % 1.0 for i in range(n)]
    return [plt.cm.hsv(h) for h in hues]


def main():
    args = parse_arguments()
    if args.ntop < 0:
        raise SystemExit("--ntop must be >= 0")

    data = np.load(args.corData)
    if "matrix" in data:
        mat = data["matrix"]
    elif "data" in data:
        mat = data["data"]
    else:
        raise SystemExit("npz file must contain 'matrix' (or 'data') and 'labels' arrays")
    labels = [str(x) for x in data["labels"]]
    if mat.ndim != 2 or mat.shape[1] != len(labels):
        raise SystemExit("matrix shape %s does not match %d labels" % (mat.shape, len(labels)))
    m = np.asarray(mat, dtype=float)

    # filter rows
    rvs = m.var(axis=1)
    keep = np.nonzero(rvs)[0]                       # drop zero-variance rows
    m = m[keep, :]
    rvs = rvs[keep]
    if args.ntop > 0 and m.shape[0] > args.ntop:
        m = m[np.argpartition(rvs, -args.ntop)[-args.ntop:], :]

    # center + scale each sample column (z-score)
    m2 = (m - np.mean(m, axis=0))
    m2 /= np.std(m2, axis=0, ddof=1)

    # SVD
    U, s, Vh = np.linalg.svd(m2, full_matrices=False)
    eigenvalues = s**2
    variance = eigenvalues / float(max(1, m2.shape[1] - 1))
    pvar = variance / variance.sum()

    # write data table (same layout as plotPCA --outFileNameData)
    if args.outFileNameData:
        with open(args.outFileNameData, "w") as of:
            of.write("#plot_pca.py --outFileNameData\n")
            of.write("Component\t%s\tEigenvalue\n" % "\t".join(labels))
            for i in range(eigenvalues.size):
                of.write("%d\t%s\t%s\n" % (i + 1, "\t".join("%s" % x for x in Vh[i, :]), eigenvalues[i]))

    # group assignment: samplesheet group column, prefix fallback per sample
    sample_groups = {}
    if args.samplesheet:
        try:
            sample_groups = load_sample_groups(args.samplesheet)
        except Exception as e:
            print("WARNING: could not read samplesheet %s (%s); using prefix groups" % (args.samplesheet, e))
            sample_groups = {}

    groups = {}
    for i, label in enumerate(labels):
        g = sample_groups.get(label)
        if not g:
            g = group_prefix(label)
            if args.samplesheet:
                print("WARNING: sample %s not in samplesheet; using prefix group '%s'" % (label, g))
        groups.setdefault(g, []).append((Vh[0, i], Vh[1, i], label))

    # plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(args.plotWidth, args.plotHeight),
                                   gridspec_kw={"height_ratios": [3, 1]})
    ax1.axhline(y=0, color="black", linestyle="dotted", zorder=1)
    ax1.axvline(x=0, color="black", linestyle="dotted", zorder=2)

    colors = categorical_colors(len(groups))
    for (g, pts), color in zip(sorted(groups.items()), colors):
        xs = [pt[0] for pt in pts]
        ys = [pt[1] for pt in pts]
        ax1.scatter(xs, ys, marker="o", color=color, s=150, label=g, zorder=3)

    # adaptive axis limits: data range + 10% padding
    def padded_limits(vals):
        vmin, vmax = min(vals), max(vals)
        span = vmax - vmin
        if span == 0:
            span = max(abs(vmax), 1.0)
        pad = 0.1 * span
        return vmin - pad, vmax + pad

    ax1.set_xlim(*padded_limits(Vh[0, :]))
    ax1.set_ylim(*padded_limits(Vh[1, :]))
    if args.plotTitle:
        ax1.set_title(args.plotTitle)
    ax1.set_xlabel("PC1 (%.1f%% of var. explained)" % (100.0 * pvar[0]))
    ax1.set_ylabel("PC2 (%.1f%% of var. explained)" % (100.0 * pvar[1]))
    ax1.legend(scatterpoints=1, loc="center left", borderaxespad=0.5,
               bbox_to_anchor=(1, 0.5), fontsize=10)

    # scree plot
    n_bars = min(len(labels), eigenvalues.size)
    ind = np.arange(n_bars)
    ax2.bar(ind, eigenvalues[:n_bars], width=0.7)
    ax2.set_ylabel("Eigenvalue")
    ax2.set_xlabel("Principal Component")
    ax2.set_xticks(ind)
    ax2.set_xticklabels(ind + 1)
    ax3 = ax2.twinx()
    ax3.plot(ind, pvar.cumsum()[:n_bars], "r-", marker="o", markersize=4)
    ax3.set_ylim([0, 1.05])
    ax3.set_ylabel("Cumulative variability")

    plt.tight_layout()
    plt.savefig(args.plotFile, bbox_inches="tight")
    plt.close()

    print("Wrote %s (PC1=%.1f%%, PC2=%.1f%%, %d groups)" % (args.plotFile, 100 * pvar[0], 100 * pvar[1], len(groups)))


if __name__ == "__main__":
    main()
