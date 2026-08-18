#!/usr/bin/env python3
"""Render manuscript Figure 4 and standalone panels from per-locus data."""

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
MISMATCH_SOURCE_DATA_PATH = FIGURE_DIR / "figure4_mismatch_source_data.tsv"
SUMMARY_PATH = BENCHMARK_DIR / "correctness_summary.tsv"
PDF_PATH = FIGURE_DIR / "figure4_candidate_characterization.pdf"
PNG_PATH = FIGURE_DIR / "figure4_candidate_characterization.png"
PANEL_DIR = FIGURE_DIR / "panels"

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
    "exact": ("Exact", BASELINE_COLOR),
    "guide_insertion_1nt": ("1-nt guide insertion", COLUMBA_COLOR),
    "guide_insertion_2nt": ("2-nt guide insertion", LIGHT_BLUE),
    "reference_deletion_1nt": ("1-nt reference deletion", ORANGE),
    "reference_deletion_2nt": ("2-nt reference deletion", YELLOW),
}
EVENT_TICK_LABELS = ["Exact", "1-nt I", "2-nt I", "1-nt D", "2-nt D"]
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


def validate_mismatch_source(rows: list[dict[str, str]]) -> dict[int, dict[str, str]]:
    summaries = {
        int(row["max_mismatches"]): row
        for row in rows
        if row["record_type"] == "summary"
    }
    if set(summaries) != {1, 2}:
        raise ValueError(f"Expected mismatch settings 1 and 2, found {sorted(summaries)}")

    for max_mismatches, summary in summaries.items():
        if int(summary["max_bulges"]) != 0 or int(summary["max_bulge_size"]) != 0:
            raise ValueError("Mismatch validation must not allow bulges")
        baseline = int(summary["baseline_valid_loci"])
        recovered = int(summary["columba_recovered_baseline_valid_loci"])
        if baseline != recovered:
            raise ValueError(f"m={max_mismatches}: {recovered}/{baseline} recovered")
        for field in (
            "baseline_missing_loci",
            "columba_only_valid_loci",
            "columba_invalid_loci",
            "baseline_command_failures",
            "columba_command_failures",
        ):
            if int(summary[field]) != 0:
                raise ValueError(f"m={max_mismatches}: nonzero {field}")
        if summary["baseline_deterministic"] != "yes" or summary["columba_deterministic"] != "yes":
            raise ValueError(f"m={max_mismatches}: nondeterministic output")

        mismatch_rows = [
            row
            for row in rows
            if row["record_type"] == "mismatch_count"
            and int(row["max_mismatches"]) == max_mismatches
        ]
        strand_rows = [
            row
            for row in rows
            if row["record_type"] == "strand"
            and int(row["max_mismatches"]) == max_mismatches
        ]
        expected_mismatch_counts = set(range(max_mismatches + 1))
        observed_mismatch_counts = {int(row["mismatch_count"]) for row in mismatch_rows}
        if observed_mismatch_counts != expected_mismatch_counts:
            raise ValueError(
                f"m={max_mismatches}: mismatch counts {sorted(observed_mismatch_counts)}"
            )
        if {row["strand"] for row in strand_rows} != {"+", "-"}:
            raise ValueError(f"m={max_mismatches}: missing strand summary")

        for grouped_rows, label in ((mismatch_rows, "mismatch"), (strand_rows, "strand")):
            grouped_baseline = sum(int(row["baseline_valid_loci"]) for row in grouped_rows)
            grouped_recovered = sum(
                int(row["columba_recovered_baseline_valid_loci"])
                for row in grouped_rows
            )
            grouped_missing = sum(int(row["baseline_missing_loci"]) for row in grouped_rows)
            if (grouped_baseline, grouped_recovered, grouped_missing) != (
                baseline,
                recovered,
                0,
            ):
                raise ValueError(
                    f"m={max_mismatches}: inconsistent {label} totals "
                    f"{grouped_baseline}/{grouped_recovered}/{grouped_missing}"
                )
    return summaries


def intended_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row["intended_valid"] == "yes"]


def event_counts(rows: list[dict[str, str]]) -> dict[str, Counter[str]]:
    valid = intended_rows(rows)
    return {
        group: Counter(
            row["canonical_event_class"]
            for row in valid
            if row["comparison_class"] == group
        )
        for group in COMPARISON_STYLE
    }


def edit_distance_counts(rows: list[dict[str, str]]) -> dict[str, Counter[int]]:
    valid = intended_rows(rows)
    return {
        group: Counter(
            int(row["canonical_edit_distance"])
            for row in valid
            if row["comparison_class"] == group
        )
        for group in COMPARISON_STYLE
    }


def pam_counts(rows: list[dict[str, str]]) -> Counter[str]:
    return Counter(row["observed_pam"] for row in rows)


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

    valid = intended_rows(rows)
    if any(row["requested_pam_match"] != "yes" for row in valid):
        raise ValueError("An intended-valid row does not have the requested PAM")
    if any(int(row["canonical_mismatches"]) != 0 for row in valid):
        raise ValueError("Canonical mismatches are nonzero despite m=0")
    if any(row["canonical_event_class"] not in EVENT_STYLE for row in valid):
        raise ValueError("Unexpected canonical event class")
    if any(int(row["canonical_edit_distance"]) not in (0, 1, 2) for row in valid):
        raise ValueError("Unexpected canonical edit distance")

    events = event_counts(rows)
    distances = edit_distance_counts(rows)
    for group in COMPARISON_STYLE:
        event_total = sum(events[group].values())
        distance_total = sum(distances[group].values())
        expected = counts[
            "shared_baseline_loci"
            if group == "shared_baseline"
            else "columba_only_valid_loci"
        ]
        if event_total != expected or distance_total != expected:
            raise ValueError(
                f"{group}: event total {event_total}, ED total {distance_total}, expected {expected}"
            )

    pams = pam_counts(rows)
    gg = pams["GG"]
    non_gg = len(rows) - gg
    if gg != counts["columba_valid_loci"]:
        raise ValueError(f"GG count {gg} != valid-locus count {counts['columba_valid_loci']}")
    if non_gg != counts["columba_non_gg_pam_records"]:
        raise ValueError(
            f"non-GG count {non_gg} != summary {counts['columba_non_gg_pam_records']}"
        )
    if counts["shared_baseline_loci"] + counts["columba_only_valid_loci"] != gg:
        raise ValueError("Shared and Columba-only counts do not sum to requested-GG loci")
    return counts


def apply_style(standalone: bool = False) -> None:
    base = 9.5 if standalone else 7.5
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": base,
            "axes.labelsize": base,
            "axes.titlesize": 11.0 if standalone else 8.5,
            "axes.titleweight": "semibold",
            "axes.linewidth": 0.7,
            "xtick.labelsize": 8.8 if standalone else 7.0,
            "ytick.labelsize": 8.8 if standalone else 7.0,
            "legend.fontsize": 8.8 if standalone else 6.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def configure_axis(axis) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#555555")
    axis.spines["bottom"].set_color("#555555")
    axis.tick_params(width=0.7, length=3)
    axis.grid(axis="y", color="#E1E1E1", linewidth=0.55, zorder=0)


def annotate_bars(axis, bars, values: list[int], fontsize: float) -> None:
    for bar, value in zip(bars, values, strict=True):
        axis.annotate(
            f"{value:,}",
            (bar.get_x() + bar.get_width() / 2, max(value, 0)),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=fontsize,
            color=DARK_TEXT,
        )


def panel_a(axis, rows: list[dict[str, str]], standalone: bool = False) -> None:
    counts = event_counts(rows)
    positions = list(range(len(EVENT_ORDER)))
    width = 0.36
    for offset, group in ((-width / 2, "shared_baseline"), (width / 2, "columba_only")):
        values = [counts[group][event] for event in EVENT_ORDER]
        label, color = COMPARISON_STYLE[group]
        bars = axis.bar(
            [position + offset for position in positions],
            values,
            width=width,
            color=color,
            edgecolor="white",
            linewidth=0.5,
            label=f"{label} (n={sum(values):,})",
            zorder=3,
        )
        annotate_bars(axis, bars, values, 7.5 if standalone else 5.8)

    axis.set_xticks(positions, EVENT_TICK_LABELS, rotation=20, ha="right")
    axis.set_ylabel("Number of loci")
    axis.set_title("Edit/event composition", loc="left", pad=15 if not standalone else 18)
    axis.legend(loc="upper right", frameon=False)
    axis.margins(y=0.12)
    axis.text(
        0.0,
        1.015,
        "Canonical oracle representation; m=0",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.0 if standalone else 6.1,
        color="#444444",
    )
    axis.text(
        0.0,
        -0.28 if standalone else -0.23,
        "I: guide insertion relative to reference; D: reference-consuming deletion",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=7.8 if standalone else 5.9,
        color="#444444",
    )


def panel_b(axis, rows: list[dict[str, str]], standalone: bool = False) -> None:
    bases = "ACGT"
    dinucleotides = [left + right for left in bases for right in bases]
    counts = pam_counts(rows)
    values = [counts[pam] for pam in dinucleotides]
    colors = [COLUMBA_COLOR if pam == "GG" else LIGHT_GRAY for pam in dinucleotides]
    bars = axis.bar(range(len(dinucleotides)), values, color=colors, width=0.78, zorder=3)
    axis.set_yscale("log")
    axis.set_xticks(range(len(dinucleotides)), dinucleotides, rotation=45, ha="right")
    axis.set_xlabel("Observed PAM dinucleotide")
    axis.set_ylabel("Final PAF records (log scale)")
    axis.set_title("Observed PAM composition", loc="left", pad=7)
    gg = counts["GG"]
    non_gg = len(rows) - gg
    axis.text(
        0.98,
        0.95,
        f"Requested GG: {gg:,}\nnon-GG: {non_gg:,}\ntotal: {len(rows):,}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.8 if standalone else 6.8,
        color=DARK_TEXT,
    )
    if standalone:
        gg_index = dinucleotides.index("GG")
        bar = bars[gg_index]
        axis.annotate(
            f"{gg:,}",
            (bar.get_x() + bar.get_width() / 2, gg),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8.2,
            color=COLUMBA_COLOR,
        )


def panel_c(
    axis,
    mismatch_rows: list[dict[str, str]],
    standalone: bool = False,
) -> None:
    summaries = {
        int(row["max_mismatches"]): row
        for row in mismatch_rows
        if row["record_type"] == "summary"
    }
    summary_m2 = summaries[2]
    summary_m1 = summaries[1]
    rows_m2 = sorted(
        (
            row
            for row in mismatch_rows
            if row["record_type"] == "mismatch_count"
            and int(row["max_mismatches"]) == 2
        ),
        key=lambda row: int(row["mismatch_count"]),
    )
    positions = [int(row["mismatch_count"]) for row in rows_m2]
    baseline_values = [int(row["baseline_valid_loci"]) for row in rows_m2]
    recovered_values = [
        int(row["columba_recovered_baseline_valid_loci"]) for row in rows_m2
    ]
    width = 0.36
    baseline_bars = axis.bar(
        [position - width / 2 for position in positions],
        baseline_values,
        width=width,
        color=BASELINE_COLOR,
        edgecolor="white",
        linewidth=0.5,
        label="Baseline-valid loci",
        zorder=3,
    )
    recovered_bars = axis.bar(
        [position + width / 2 for position in positions],
        recovered_values,
        width=width,
        color=COLUMBA_COLOR,
        edgecolor="white",
        linewidth=0.5,
        label="Recovered by Columba",
        zorder=3,
    )
    annotate_bars(axis, baseline_bars, baseline_values, 7.5 if standalone else 5.8)
    annotate_bars(axis, recovered_bars, recovered_values, 7.5 if standalone else 5.8)

    axis.set_yscale("log")
    axis.set_ylim(1, 2500)
    axis.set_xticks(positions, [str(position) for position in positions])
    axis.set_xlabel("Mismatch count")
    axis.set_ylabel("Number of loci (log scale)")
    axis.set_title("Real-reference mismatch recovery", loc="left", pad=15 if not standalone else 18)
    axis.legend(loc="upper left", frameon=False)
    axis.text(
        0.0,
        1.015,
        "chr22, 10 guides; m=2, b=0, z=0",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.0 if standalone else 6.1,
        color="#444444",
    )
    axis.text(
        0.03,
        0.62,
        f"{int(summary_m2['columba_recovered_baseline_valid_loci']):,}/"
        f"{int(summary_m2['baseline_valid_loci']):,} recovered\n"
        f"{int(summary_m2['baseline_missing_loci']):,} missing; "
        f"{int(summary_m2['columba_invalid_loci']):,} invalid\n"
        f"m=1: {int(summary_m1['columba_recovered_baseline_valid_loci']):,}/"
        f"{int(summary_m1['baseline_valid_loci']):,} recovered",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.0 if standalone else 6.1,
        color=DARK_TEXT,
    )


def draw_box(
    axis,
    xy,
    width,
    height,
    text,
    facecolor,
    edgecolor,
    fontsize: float,
    text_color=DARK_TEXT,
):
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
        fontsize=fontsize,
        color=text_color,
        fontweight="semibold",
    )


def draw_arrow(axis, start, end, scale: float = 8) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        transform=axis.transAxes,
        arrowstyle="-|>",
        mutation_scale=scale,
        linewidth=0.8,
        color="#666666",
        connectionstyle="arc3,rad=0.0",
    )
    axis.add_patch(arrow)


def panel_d(
    axis,
    rows: list[dict[str, str]],
    counts: dict[str, int],
    standalone: bool = False,
) -> None:
    axis.set_axis_off()
    axis.set_title("Post-output classification", loc="left", pad=8)
    final_count = len(rows)
    oracle_valid = sum(row["oracle_valid"] == "yes" for row in rows)
    gg_valid = sum(row["intended_valid"] == "yes" for row in rows)
    shared = counts["shared_baseline_loci"]
    columba_only = counts["columba_only_valid_loci"]
    non_gg = counts["columba_non_gg_pam_records"]
    invalid = counts["columba_invalid_records"]
    box_font = 9.5 if standalone else 7.0

    axis.text(
        0.0,
        .965,
        "Classification of final output; not a candidate-generation or filtering funnel",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.4 if standalone else 6.0,
        color="#444444",
    )
    draw_box(
        axis,
        (0.34, 0.79),
        0.32,
        0.12,
        f"Final Columba PAF records\n{final_count:,}",
        "#F2F2F2",
        "#777777",
        box_font,
    )
    draw_box(
        axis,
        (0.31, 0.58),
        0.38,
        0.12,
        f"Independently alignment-valid\n{oracle_valid:,}",
        "#E8F3F8",
        COLUMBA_COLOR,
        box_font,
    )
    draw_box(
        axis,
        (0.34, 0.37),
        0.32,
        0.12,
        f"Requested-GG/oracle-valid loci\n{gg_valid:,}",
        "#D9EEF8",
        COLUMBA_COLOR,
        box_font,
    )
    draw_arrow(axis, (0.50, 0.79), (0.50, 0.70), 10 if standalone else 8)
    draw_arrow(axis, (0.50, 0.58), (0.50, 0.49), 10 if standalone else 8)
    axis.text(
        0.76,
        0.63,
        f"{invalid:,} invalid",
        transform=axis.transAxes,
        ha="left",
        fontsize=8.2 if standalone else 6.2,
        color="#555555",
    )
    axis.text(
        0.76,
        0.42,
        f"{non_gg:,} non-GG",
        transform=axis.transAxes,
        ha="left",
        fontsize=8.2 if standalone else 6.2,
        color="#555555",
    )

    draw_box(
        axis,
        (0.03, 0.07),
        0.39,
        0.16,
        f"Shared baseline loci\n{shared:,}",
        "#E1E1E1",
        BASELINE_COLOR,
        box_font,
    )
    draw_box(
        axis,
        (0.52, 0.07),
        0.45,
        0.16,
        f"Additional oracle-valid\nColumba candidate loci\n{columba_only:,}",
        "#D9EEF8",
        COLUMBA_COLOR,
        box_font,
    )
    draw_arrow(axis, (0.44, 0.37), (0.25, 0.23), 10 if standalone else 8)
    draw_arrow(axis, (0.56, 0.37), (0.74, 0.23), 10 if standalone else 8)


def save_figure(figure, pdf_path: Path, png_path: Path, subject: str) -> None:
    pdf_metadata = {
        "Title": subject,
        "Subject": subject,
        "Creator": "make_figure4.py",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(pdf_path, format="pdf", metadata=pdf_metadata, facecolor="white")
    figure.savefig(
        png_path,
        format="png",
        dpi=320,
        metadata={"Software": "make_figure4.py"},
        facecolor="white",
    )
    plt.close(figure)


def render_combined(
    rows: list[dict[str, str]],
    counts: dict[str, int],
    mismatch_rows: list[dict[str, str]],
) -> None:
    apply_style(standalone=False)
    figure, axes = plt.subplots(2, 2, figsize=(8.3, 6.2))
    figure.subplots_adjust(
        left=0.083,
        right=0.985,
        bottom=0.105,
        top=0.94,
        wspace=0.33,
        hspace=0.49,
    )
    axis_a, axis_b, axis_c, axis_d = axes.flat
    for axis in (axis_a, axis_b, axis_c):
        configure_axis(axis)

    panel_a(axis_a, rows)
    panel_b(axis_b, rows)
    panel_c(axis_c, mismatch_rows)
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

    save_figure(
        figure,
        PDF_PATH,
        PNG_PATH,
        "Figure 4: Characterization of Columba-reported candidate loci",
    )


def render_standalone(
    rows: list[dict[str, str]],
    counts: dict[str, int],
    mismatch_rows: list[dict[str, str]],
    panel_name: str,
    render_panel,
    figsize: tuple[float, float],
    margins: dict[str, float],
) -> tuple[Path, Path]:
    apply_style(standalone=True)
    figure, axis = plt.subplots(figsize=figsize)
    figure.subplots_adjust(**margins)
    if panel_name != "panelD_post_output_classification":
        configure_axis(axis)
    if panel_name == "panelD_post_output_classification":
        render_panel(axis, rows, counts, standalone=True)
    elif panel_name == "panelC_mismatch_validation":
        render_panel(axis, mismatch_rows, standalone=True)
    else:
        render_panel(axis, rows, standalone=True)
    pdf_path = PANEL_DIR / f"{panel_name}.pdf"
    png_path = PANEL_DIR / f"{panel_name}.png"
    save_figure(figure, pdf_path, png_path, panel_name.replace("_", " "))
    return pdf_path, png_path


def render_all(
    rows: list[dict[str, str]],
    counts: dict[str, int],
    mismatch_rows: list[dict[str, str]],
) -> list[Path]:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    render_combined(rows, counts, mismatch_rows)
    outputs = [PDF_PATH, PNG_PATH]
    standalone_specs = [
        (
            "panelA_edit_event_composition",
            panel_a,
            (6.8, 4.5),
            {"left": 0.12, "right": 0.98, "bottom": 0.25, "top": 0.87},
        ),
        (
            "panelB_pam_distribution",
            panel_b,
            (6.8, 4.2),
            {"left": 0.12, "right": 0.98, "bottom": 0.18, "top": 0.90},
        ),
        (
            "panelC_mismatch_validation",
            panel_c,
            (6.2, 4.3),
            {"left": 0.13, "right": 0.98, "bottom": 0.14, "top": 0.86},
        ),
        (
            "panelD_post_output_classification",
            panel_d,
            (6.6, 4.8),
            {"left": 0.05, "right": 0.97, "bottom": 0.06, "top": 0.89},
        ),
    ]
    for spec in standalone_specs:
        outputs.extend(render_standalone(rows, counts, mismatch_rows, *spec))
    return outputs


def main() -> None:
    rows = read_tsv(SOURCE_DATA_PATH)
    mismatch_rows = read_tsv(MISMATCH_SOURCE_DATA_PATH)
    counts = validate_source(rows)
    validate_mismatch_source(mismatch_rows)
    outputs = render_all(rows, counts, mismatch_rows)
    for path in outputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")
        print(f"wrote {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
