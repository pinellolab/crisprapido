#!/usr/bin/env python3
"""Render manuscript Figure 4 from compact per-locus source data."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


FIGURE_DIR = Path(__file__).resolve().parent
REPO_ROOT = FIGURE_DIR.parents[2]
BENCHMARK_DIR = (
    REPO_ROOT
    / "benchmark"
    / "performance"
    / "chm13v2_whole_genome_20_guides"
)
SOURCE_DATA_PATH = FIGURE_DIR / "figure4_source_data.tsv"
SUMMARY_PATH = BENCHMARK_DIR / "correctness_summary.tsv"
PDF_PATH = FIGURE_DIR / "figure4_candidate_characterization.pdf"
PNG_PATH = FIGURE_DIR / "figure4_candidate_characterization.png"

BASELINE_COLOR = "#4D4D4D"
COLUMBA_COLOR = "#0072B2"
LIGHT_BLUE = "#56B4E9"
ORANGE = "#D55E00"
YELLOW = "#E69F00"
LIGHT_GRAY = "#C7C7C7"
DARK_TEXT = "#262626"

EVENT_ORDER = [
    "exact",
    "guide_insertion_1nt",
    "guide_insertion_2nt",
    "reference_deletion_1nt",
    "reference_deletion_2nt",
]
EVENT_STYLE = {
    "exact": ("Exact (ED0)", BASELINE_COLOR),
    "guide_insertion_1nt": ("1-nt I (ED1)", COLUMBA_COLOR),
    "guide_insertion_2nt": ("2-nt I (ED2)", LIGHT_BLUE),
    "reference_deletion_1nt": ("1-nt D (ED1)", ORANGE),
    "reference_deletion_2nt": ("2-nt D (ED2)", YELLOW),
}
COMPARISON_STYLE = {
    "shared_baseline": ("Shared with baseline", BASELINE_COLOR),
    "columba_only": ("Columba-only", COLUMBA_COLOR),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No rows in {path}")
    return rows


def validate_source(rows: list[dict[str, str]]) -> dict[str, int]:
    summary = read_tsv(SUMMARY_PATH)[0]
    counts = {
        "columba_raw_records": len(rows),
        "columba_valid_loci": sum(row["intended_valid"] == "yes" for row in rows),
        "shared_baseline_loci": sum(
            row["comparison_class"] == "shared_baseline" for row in rows
        ),
        "columba_only_valid_loci": sum(
            row["comparison_class"] == "columba_only" for row in rows
        ),
        "columba_invalid_records": sum(row["oracle_valid"] == "no" for row in rows),
        "columba_non_gg_pam_records": sum(
            row["oracle_valid"] == "yes" and row["requested_pam_match"] == "no"
            for row in rows
        ),
    }
    counts["baseline_valid_loci"] = int(summary["baseline_valid_loci"])
    counts["baseline_missing_from_columba"] = (
        counts["baseline_valid_loci"] - counts["shared_baseline_loci"]
    )
    for field, observed in counts.items():
        expected = int(summary[field])
        if observed != expected:
            raise ValueError(f"{field}: source data {observed}, summary {expected}")

    intended = [row for row in rows if row["intended_valid"] == "yes"]
    if any(row["requested_pam_match"] != "yes" for row in intended):
        raise ValueError("An intended-valid row does not have the requested PAM")
    if any(int(row["canonical_mismatches"]) != 0 for row in intended):
        raise ValueError("Canonical mismatches are nonzero despite m=0")
    if any(row["canonical_event_class"] not in EVENT_STYLE for row in intended):
        raise ValueError("Unexpected canonical event class")
    if any(not (0.0 <= float(row["cfd_score"]) <= 1.0) for row in rows):
        raise ValueError("CFD score outside [0, 1]")
    return counts


def configure_axis(axis) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#555555")
    axis.spines["bottom"].set_color("#555555")
    axis.tick_params(width=0.7, length=3)
    axis.grid(axis="y", color="#E1E1E1", linewidth=0.55, zorder=0)


def panel_a(axis, rows: list[dict[str, str]]) -> None:
    classes = ["shared_baseline", "columba_only"]
    labels = ["Shared", "Columba-only"]
    positions = [0, 1]
    bottoms = [0.0, 0.0]
    totals = [sum(row["comparison_class"] == group for row in rows) for group in classes]

    for event in EVENT_ORDER:
        values = []
        for group, total in zip(classes, totals, strict=True):
            count = sum(
                row["comparison_class"] == group
                and row["canonical_event_class"] == event
                for row in rows
            )
            values.append(100.0 * count / total)
        bars = axis.bar(
            positions,
            values,
            bottom=bottoms,
            width=0.62,
            color=EVENT_STYLE[event][1],
            edgecolor="white",
            linewidth=0.6,
            label=EVENT_STYLE[event][0],
            zorder=3,
        )
        for index, (bar, value) in enumerate(zip(bars, values, strict=True)):
            if value >= 7.0:
                axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    bottoms[index] + value / 2,
                    f"{value:.0f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=6.4,
                    fontweight="semibold",
                )
            bottoms[index] += value

    axis.set_xticks(positions, [f"{label}\n(n={total:,})" for label, total in zip(labels, totals, strict=True)])
    axis.set_ylim(0, 100)
    axis.set_ylabel("Composition of oracle-valid GG loci (%)")
    axis.set_title("Edit/event composition", loc="left", pad=16)
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.17),
        ncol=3,
        frameon=False,
        fontsize=6.1,
        handlelength=1.2,
        columnspacing=1.0,
    )
    axis.text(
        0.0,
        1.015,
        "Canonical oracle representation; m=0 (no substitution distribution)",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.1,
        color="#444444",
    )


def panel_b(axis, rows: list[dict[str, str]]) -> None:
    bases = "ACGT"
    dinucleotides = [left + right for left in bases for right in bases]
    counts = Counter(row["observed_pam"] for row in rows)
    values = [counts[pam] for pam in dinucleotides]
    colors = [COLUMBA_COLOR if pam == "GG" else LIGHT_GRAY for pam in dinucleotides]
    axis.bar(range(len(dinucleotides)), values, color=colors, width=0.78, zorder=3)
    axis.set_yscale("log")
    axis.set_xticks(range(len(dinucleotides)), dinucleotides, rotation=45, ha="right")
    axis.set_ylabel("Final PAF records (log scale)")
    axis.set_title("Observed PAM composition", loc="left", pad=7)
    gg = counts["GG"]
    non_gg = len(rows) - gg
    axis.text(
        0.98,
        0.95,
        f"GG: {gg:,}\nnon-GG: {non_gg:,}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=6.8,
        color=DARK_TEXT,
    )


def ecdf(values: list[float]) -> tuple[list[float], list[float]]:
    ordered = sorted(values)
    cumulative = [100.0 * index / len(ordered) for index in range(1, len(ordered) + 1)]
    return ordered, cumulative


def panel_c(axis, rows: list[dict[str, str]]) -> None:
    for group in ("shared_baseline", "columba_only"):
        values = [
            float(row["cfd_score"])
            for row in rows
            if row["comparison_class"] == group
        ]
        x_values, y_values = ecdf(values)
        label, color = COMPARISON_STYLE[group]
        axis.step(
            x_values,
            y_values,
            where="post",
            color=color,
            linewidth=1.7,
            label=f"{label} (n={len(values):,})",
        )
    axis.set_xlim(0, 1.0)
    axis.set_ylim(0, 100)
    axis.set_xlabel("Reported CFD score")
    axis.set_ylabel("Cumulative loci (%)")
    axis.set_title("CFD scores for oracle-valid GG loci", loc="left", pad=16)
    axis.legend(loc="lower right", frameon=False, fontsize=6.7)
    axis.text(
        0.0,
        1.015,
        "Non-GG records excluded; scores retain PAF precision",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.1,
        color="#444444",
    )


def draw_box(axis, xy, width, height, text, facecolor, edgecolor, text_color=DARK_TEXT):
    box = Rectangle(
        xy,
        width,
        height,
        transform=axis.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=1.1,
    )
    axis.add_patch(box)
    axis.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=7.0,
        color=text_color,
        fontweight="semibold",
    )


def draw_arrow(axis, start, end) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        transform=axis.transAxes,
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=0.8,
        color="#666666",
        connectionstyle="arc3,rad=0.0",
    )
    axis.add_patch(arrow)


def panel_d(axis, rows: list[dict[str, str]], counts: dict[str, int]) -> None:
    axis.set_axis_off()
    axis.set_title(
        "Post-output classification\n(not a candidate/filtering funnel)",
        loc="left",
        pad=4,
        fontsize=8.0,
    )
    final_count = len(rows)
    oracle_valid = sum(row["oracle_valid"] == "yes" for row in rows)
    gg_valid = sum(row["intended_valid"] == "yes" for row in rows)
    shared = counts["shared_baseline_loci"]
    columba_only = counts["columba_only_valid_loci"]
    non_gg = counts["columba_non_gg_pam_records"]
    invalid = counts["columba_invalid_records"]

    draw_box(axis, (0.34, 0.79), 0.32, 0.12, f"Final PAF\n{final_count:,}", "#F2F2F2", "#777777")
    draw_box(axis, (0.31, 0.58), 0.38, 0.12, f"Oracle alignment-valid\n{oracle_valid:,}", "#E8F3F8", COLUMBA_COLOR)
    draw_box(axis, (0.36, 0.37), 0.28, 0.12, f"Requested GG\n{gg_valid:,}", "#D9EEF8", COLUMBA_COLOR)
    draw_arrow(axis, (0.50, 0.79), (0.50, 0.70))
    draw_arrow(axis, (0.50, 0.58), (0.50, 0.49))
    axis.text(0.71, 0.63, f"{invalid:,} invalid", transform=axis.transAxes, ha="left", fontsize=6.2, color="#555555")
    axis.text(0.67, 0.42, f"{non_gg:,} non-GG", transform=axis.transAxes, ha="left", fontsize=6.2, color="#555555")

    draw_box(axis, (0.04, 0.08), 0.38, 0.15, f"Shared baseline loci\n{shared:,}", "#E1E1E1", BASELINE_COLOR)
    draw_box(axis, (0.53, 0.08), 0.43, 0.15, f"Additional Columba-only\noracle-valid loci\n{columba_only:,}", "#D9EEF8", COLUMBA_COLOR)
    draw_arrow(axis, (0.44, 0.37), (0.25, 0.23))
    draw_arrow(axis, (0.56, 0.37), (0.74, 0.23))


def render(rows: list[dict[str, str]], counts: dict[str, int]) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.5,
            "axes.titleweight": "semibold",
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure, axes = plt.subplots(2, 2, figsize=(8.3, 6.2))
    figure.subplots_adjust(left=0.083, right=0.985, bottom=0.11, top=0.94, wspace=0.33, hspace=0.47)
    axis_a, axis_b, axis_c, axis_d = axes.flat
    for axis in (axis_a, axis_b, axis_c):
        configure_axis(axis)

    panel_a(axis_a, rows)
    panel_b(axis_b, rows)
    panel_c(axis_c, rows)
    panel_d(axis_d, rows, counts)

    for label, axis in zip(("A", "B", "C", "D"), axes.flat, strict=True):
        axis.text(
            -0.15,
            1.13,
            label,
            transform=axis.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
            ha="left",
        )

    pdf_metadata = {
        "Title": "Figure 4: Characterization of Columba-reported candidate loci",
        "Subject": "Edit events, PAMs, CFD scores, and post-output classification",
        "Creator": "make_figure4.py",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(PDF_PATH, format="pdf", metadata=pdf_metadata, facecolor="white")
    figure.savefig(
        PNG_PATH,
        format="png",
        dpi=320,
        metadata={"Software": "make_figure4.py"},
        facecolor="white",
    )
    plt.close(figure)


def main() -> None:
    rows = read_tsv(SOURCE_DATA_PATH)
    counts = validate_source(rows)
    render(rows, counts)
    for path in (PDF_PATH, PNG_PATH):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")
        print(f"wrote {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
