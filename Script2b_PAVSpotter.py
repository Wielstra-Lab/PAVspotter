#!/usr/bin/env python3
# Genetic Data Analysis and Visualization
# =========================================================================
# Processes genetic data across species, analyzing gene presence/absence. 
# Supports customizable parameters for targeted analysis, including
# data summarization, visualization, and cross-correlation analysis. 
# Generates detailed plots and exports results to CSV.
#
# Usage:
# - Uncomment certain lines for standalone use in MATLAB.
# - Adjust parameters in the 'Parameters' section as needed in standalone 
#   use.
# - Ensure correct working directory.
#
# Outputs:
# - Visual plots in JPEG format (other formats supported, see Parameters.
# - CSV summary of gene analysis and cross-correlations.
#
# Requirements: Python 3.9+, input data in specified folder structure.
#
# Author: Chris van der Ploeg & Manon de Visser | Date: 2025-12-21 | 
# Contact: vanderploeg.cj@gmail.com
# =========================================================================

import os
import re
import sys
import argparse
from datetime import datetime
from typing import List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Utilities
# -----------------------------

def safe_read_lines(path: str) -> List[str]:
    """Read a text file as a list of stripped lines."""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        return [ln.strip() for ln in f.readlines() if ln.strip()]


def read_individuals(path: str, categories: List[str]) -> List[str]:
    """Read individuals.txt and keep entries containing any of the base categories."""
    lines = safe_read_lines(path)
    # Keep lines that contain any of the categories
    kept = []
    for ln in lines:
        if any(c in ln for c in categories):
            kept.append(ln)
    return kept


def read_gene_file(path: str) -> pd.DataFrame:
    """Read a gene file expected to have 4 columns: gene_name, position(index), reads, sample_name.
    The original MATLAB uses readcell(FileType='text'); this parser tries common separators.
    """
    # Try several separators; fall back to whitespace
    seps = [",", "\t", " "]
    for sep in seps:
        try:
            df = pd.read_csv(path, sep=sep, header=None, engine="python")
            if df.shape[1] >= 4:
                break
        except Exception:
            continue
    else:
        # Final attempt: any whitespace
        df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    # Ensure we have at least 4 columns
    if df.shape[1] < 4:
        raise ValueError(f"Gene file '{path}' must have at least 4 columns; got {df.shape[1]}")
    df = df.iloc[:, :4].copy()
    df.columns = ["gene_name", "position", "reads", "sample_name"]
    # Coerce types
    df["gene_name"] = df["gene_name"].astype(str)
    # Position is often integer 1..N
    df["position"] = pd.to_numeric(df["position"], errors="coerce").fillna(0).astype(int)
    df["reads"] = pd.to_numeric(df["reads"], errors="coerce").fillna(0.0).astype(float)
    df["sample_name"] = df["sample_name"].astype(str)
    return df


def xcorr_full(a: np.ndarray, b: np.ndarray, normalized: bool) -> np.ndarray:
    """Full cross-correlation (like MATLAB xcorr), optionally normalized.
    Returns a vector of length 2*n-1 with lags from -(n-1)..(n-1)."""
    if a.ndim != 1 or b.ndim != 1:
        a = a.ravel()
        b = b.ravel()
    n = len(a)
    if len(b) != n:
        raise ValueError("xcorr_full requires equal-length sequences")
    corr = np.correlate(a, b, mode="full")
    if normalized:
        na = np.linalg.norm(a)
        nb = np.linalg.norm(b)
        if na > 0 and nb > 0:
            corr = corr / (na * nb)
        # else leave unnormalized to avoid divide-by-zero
    return corr


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


# -----------------------------
# Core translation
# -----------------------------

def process_species(
    sidx: int,
    sidx_tot: int,
    species_path: str,
    base_categories: List[str],
    ctrl_category: str,
    common_identifier: str,
    sum_for_each_subspecies: bool,
    contig_width: int,
    reads_threshold: float,
    plot_figures: bool,
    figures_root: str,
) -> Tuple[List[List[str]], List[float]]:
    """Process one species folder, returning total_results rows and control xcorr entries."""
    total_results: List[List[str]] = []
    ctrl_xcorr_entries: List[float] = []

    individuals_path = os.path.join(species_path, "individuals.txt")
    if not os.path.isfile(individuals_path):
        raise FileNotFoundError(f"Missing individuals.txt in {species_path}")

    categories_tot = read_individuals(individuals_path, base_categories)
    if not categories_tot:
        raise ValueError(
            f"No individuals in '{individuals_path}' matched the selected categories {base_categories}"
        )

    # List gene files containing the common identifier
    all_files = sorted(os.listdir(species_path))
    gene_files = [f for f in all_files if common_identifier in f]
    

    if not gene_files:
        # No files -> nothing to do
        return total_results, ctrl_xcorr_entries

    # Prepare plotting directories
    fig_dir = None
    if plot_figures:
        fig_dir = figures_root

    # Iterate genes
    for gidx, gf in enumerate(gene_files, start=1):
        print(f"working on species {sidx} of {sidx_tot}")
        print(f"working on gene {gidx} of {len(gene_files)}")

        gf_path = os.path.join(species_path, gf)
        df = read_gene_file(gf_path)

        # Filter rows to those whose sample_name contains any entry from categories_tot
        available = []
        for ct in categories_tot:
            ct_mask = df["sample_name"].str.contains(re.escape(ct), regex=True)
            available.append(bool(ct_mask.any()))
        # Row-level mask: any of the categories_tot
        row_mask = pd.Series(False, index=df.index)
        for ct in categories_tot:
            row_mask = row_mask | df["sample_name"].str.contains(re.escape(ct), regex=True)
        df_f = df[row_mask].copy()
        if df_f.empty:
            # Nothing matched; skip this gene
            continue

        # Determine gene_length (per-sample length)
        total_rows = len(df_f)
        present_count = sum(1 for a in available if a)
        if present_count == 0:
            continue
        gene_length = total_rows // present_count
        # Sort df_f by categories_tot order and position ascending
        slices = []
        for ct in categories_tot:
            sli = df_f[df_f["sample_name"].str.contains(re.escape(ct), regex=True)].copy()
            if not sli.empty:
                sli.sort_values(by=["position"], inplace=True)
                # Trim/extend to gene_length
                if len(sli) > gene_length:
                    sli = sli.iloc[:gene_length]
                elif len(sli) < gene_length:
                    missing = pd.DataFrame({
                        "gene_name": [sli.iloc[0]["gene_name"]] * (gene_length - len(sli)) if len(sli) else [df_f.iloc[0]["gene_name"]] * (gene_length - len(sli)),
                        "position": list(range(1, gene_length + 1))[-(gene_length - len(sli)):],
                        "reads": [0.0] * (gene_length - len(sli)),
                        "sample_name": [ct] * (gene_length - len(sli)),
                    })
                    sli = pd.concat([sli, missing], ignore_index=True)
            else:
                # Not available; fabricate zero reads slice of length gene_length
                sli = pd.DataFrame({
                    "gene_name": [df_f.iloc[0]["gene_name"]] * gene_length,
                    "position": list(range(1, gene_length + 1)),
                    "reads": [0.0] * gene_length,
                    "sample_name": [ct] * gene_length,
                })
            slices.append(sli)
        df_ordered = pd.concat(slices, ignore_index=True)

        gene_name = str(df_ordered.iloc[0]["gene_name"]).replace("_", "_")
        n_samples = len(categories_tot)
        gene_data = np.zeros((gene_length, n_samples), dtype=float)
        coverage = []

        for i, ct in enumerate(categories_tot):
            sli = df_ordered.iloc[i * gene_length:(i + 1) * gene_length]
            seq = sli["reads"].to_numpy(dtype=float)
            # Apply contig-width filtering (zero-out short positive runs)
            seq_filtered = seq.copy()
            contig_counter = 0
            contig_start = 0
            for j in range(gene_length):
                if contig_counter == 0 and seq[j] > 0:
                    contig_counter = 1
                    contig_start = j
                elif contig_counter > 0 and seq[j] == 0:
                    if contig_counter < contig_width:
                        seq_filtered[contig_start:j+1] = 0.0
                    contig_counter = 0
                else:
                    seq_filtered[j] = seq[j]
                    if seq[j] > 0:
                        contig_counter += 1
                    if j == gene_length - 1 and contig_counter > 0:
                        if contig_counter < contig_width:
                            seq_filtered[contig_start:j+1] = 0.0
                        contig_counter = 0
            if np.linalg.norm(seq_filtered, ord=np.inf) < reads_threshold:
                seq_filtered[:] = 0.0
            gene_data[:, i] = seq_filtered
            coverage.append(float(np.max(seq_filtered)))

        # Merge sequences across subspecies
        if sum_for_each_subspecies:
            gene_data_tot = np.zeros((gene_length, len(base_categories)), dtype=float)
            for i, ct in enumerate(categories_tot):
                for q, base in enumerate(base_categories):
                    if base in ct:
                        gene_data_tot[:, q] += gene_data[:, i]
            # Normalize per base category
            for q in range(len(base_categories)):
                vmax = np.max(gene_data_tot[:, q])
                if vmax != 0:
                    gene_data_tot[:, q] = gene_data_tot[:, q] / vmax
            categories_n = base_categories
            gene_data_use = gene_data_tot
        else:
            gene_data_tot = gene_data.copy()
            for i in range(n_samples):
                vmax = np.max(gene_data_tot[:, i])
                if vmax != 0:
                    gene_data_tot[:, i] = gene_data_tot[:, i] / vmax
            categories_n = categories_tot
            gene_data_use = gene_data_tot

        # Control sequences
        ctrl_columns = []
        for i, ct in enumerate(categories_tot):
            if ctrl_category in ct:
                ctrl_columns.append(gene_data[:, i])
        # Controls xcorr
        for c1 in range(len(ctrl_columns)):
            for c2 in range(c1 + 1, len(ctrl_columns)):
                a = ctrl_columns[c1]
                b = ctrl_columns[c2]
                const_a = float(np.max(a)) == float(np.min(a))
                const_b = float(np.max(b)) == float(np.min(b))
                corr = xcorr_full(a, b, normalized=not (const_a and const_b))
                center = corr[len(corr)//2]
                ctrl_xcorr_entries.append(float(center))

        # Plotting
        if plot_figures:
            # Depth subplots
            fig1, axs = plt.subplots(len(categories_n), 1, figsize=(6, 8), constrained_layout=True)
            if len(categories_n) == 1:
                axs = [axs]
            for idx, cat in enumerate(categories_n):
                axs[idx].plot(gene_data_use[:, idx], linewidth=1.5)
                axs[idx].set_ylabel(f"Reads {cat}")
                if idx == len(categories_n) - 1:
                    axs[idx].set_xlabel("Target position [-]")
                if idx == 0:
                    axs[idx].set_title(gene_name)
                axs[idx].set_xlim(0, gene_length - 1)
            # X-corr figure
            fig2, ax2 = plt.subplots(figsize=(6, 6))
            ax2.set_title(gene_name)

        # Cross-correlation across categories
        plotted_pairs: List[Tuple[int, int]] = []
        for i in range(len(categories_n)):
            for j in range(len(categories_n)):
                a = gene_data_use[:, i]
                b = gene_data_use[:, j]
                const_a = float(np.max(a)) == float(np.min(a))
                const_b = float(np.max(b)) == float(np.min(b))
                corr = xcorr_full(a, b, normalized=not (const_a and const_b))
                center = corr[len(corr)//2]
                if i != j:
                    pair_sorted = tuple(sorted((i, j)))
                    if pair_sorted not in plotted_pairs:
                        plotted_pairs.append(pair_sorted)
                        if plot_figures:
                            ax2.plot(corr, linewidth=2, label=f"{categories_n[i]} -> {categories_n[j]} xcorr={center*100:.1f}%")
                        total_results.append([
                            os.path.basename(species_path),
                            gene_name,
                            categories_n[i],
                            categories_n[j],
                            f"{center:.6f}",
                            f"{max(coverage):.6f}",
                        ])

        if plot_figures:
            ax2.axvline(x=gene_length - 1, color="k", linestyle="--", linewidth=2)
            ax2.set_xlabel("Shift in target position [-]")
            ax2.set_ylabel("Normalized cross-correlation [-]")
            ax2.set_xlim(0, 2 * gene_length - 2)
            ax2.set_ylim(0, 1.5)
            ax2.legend(loc="best")
            base = df_ordered.iloc[0]["gene_name"]
            base = re.sub(r"\.Tpyg$", "", str(base))
            fig1_path = os.path.join(figures_root, f"{base}_subplots.jpg")
            fig2_path = os.path.join(figures_root, f"{base}_xcorrs.jpg")
            fig1.savefig(fig1_path, dpi=300)
            fig2.savefig(fig2_path, dpi=300)
            plt.close(fig1)
            plt.close(fig2)

    return total_results, ctrl_xcorr_entries


def run(
    root_dir: str,
    filetype: str,
    save_file_name: str,
    sum_for_each_subspecies: bool,
    plot_figures: bool,
    contig_width: int,
    reads_threshold: float,
    categories: List[str],
    ctrl_category: str,
    common_identifier: str,
):
    figures_dir = None
    if plot_figures:
        figures_dir = f"figures_{datetime.now().date()}"
        ensure_dir(figures_dir)

    species_folders = [
        d for d in sorted(os.listdir(root_dir))
        if os.path.isdir(os.path.join(root_dir, d)) and "." not in d
    ]


    all_total_results: List[List[str]] = []
    all_ctrl_xcorr: List[float] = []

    for sidx_p, species in enumerate(species_folders,start=1):
        species_path = os.path.join(root_dir, species)
        tr, ctrls = process_species(
            sidx = sidx_p,
            sidx_tot = len(species_folders),
            species_path=species_path,
            base_categories=categories,
            ctrl_category=ctrl_category,
            common_identifier=common_identifier,
            sum_for_each_subspecies=sum_for_each_subspecies,
            contig_width=contig_width,
            reads_threshold=reads_threshold,
            plot_figures=plot_figures,
            figures_root=figures_dir if figures_dir else root_dir,
        )
        all_total_results.extend(tr)
        all_ctrl_xcorr.extend(ctrls)

    if all_total_results:
        df_res = pd.DataFrame(
            all_total_results,
            columns=["species", "gene", "from_category", "to_category", "xcorr_center", "max_coverage"],
        )
        df_res.to_csv(f"{save_file_name}.csv", index=False, header=False)
    else:
        pd.DataFrame(columns=["species", "gene", "from_category", "to_category", "xcorr_center", "max_coverage"]).to_csv(
            f"{save_file_name}.csv", index=False, header=False
        )

    if all_ctrl_xcorr:
        df_ctrl = pd.DataFrame({"xcorr_center": all_ctrl_xcorr})
        df_ctrl.to_csv(f"{save_file_name}_ctrl_xcorr.csv", index=False, header=False)
        if plot_figures:
            ctrl_plot_data = (1.0 - np.array(all_ctrl_xcorr)) * 100.0
            ctrl_plot_data = np.where(ctrl_plot_data == 100.0, np.nan, ctrl_plot_data)
            plt.figure(figsize=(6, 6))
            plt.hist(ctrl_plot_data[~np.isnan(ctrl_plot_data)], bins=200)
            plt.savefig(os.path.join(figures_dir, "histogram_ctrl.jpg"), dpi=300)
            plt.close()
    else:
        pd.DataFrame({"xcorr_center": []}).to_csv(f"{save_file_name}_ctrl_xcorr.csv", index=False, header=False)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pav_spotter",
        description=(
            "Processes genetic data across species, analyzing gene presence/absence, "
            "summarization, visualization, and cross-correlation (Python conversion of MATLAB Script2b_PAVSpotter)."
        ),
    )
    p.add_argument("--root", default=".", help="Root directory containing species folders")
    p.add_argument("--filetype", default="jpg", help="Output figure file type (default: jpg)")
    p.add_argument("--save-file-name", default="all_data", help="Base name for CSV outputs")
    p.add_argument("--sum-for-each-subspecies", action="store_true", help="Sum sequences per base category across subspecies")
    p.add_argument("--no-plot", dest="plot_figures", action="store_false", help="Disable plotting")
    p.add_argument("--contig-width", type=int, default=100, help="Minimum contig width to retain positive reads")
    p.add_argument("--reads-threshold", type=float, default=10.0, help="Infinity-norm threshold below which a sequence is zeroed")
    p.add_argument(
        "--categories",
        nargs="+",
        default=["ST", "FT", "hatch"],
        help="Base categories to select (e.g., ST FT hatch)",
    )
    p.add_argument("--ctrl-category", default="hatch", help="Control category substring")
    p.add_argument("--common-identifier", default="DN", help="Filename substring to select gene files")
    p.set_defaults(plot_figures=True)
    return p


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(argv)

    run(
        root_dir=args.root,
        filetype=args.filetype,
        save_file_name=args.save_file_name,
        sum_for_each_subspecies=args.sum_for_each_subspecies,
        plot_figures=args.plot_figures,
        contig_width=args.contig_width,
        reads_threshold=args.reads_threshold,
        categories=args.categories,
        ctrl_category=args.ctrl_category,
        common_identifier=args.common_identifier,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
