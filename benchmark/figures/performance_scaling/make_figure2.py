#!/usr/bin/env python3
"""Build manuscript Figure 2 from finalized benchmark summary TSVs."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullFormatter


FIGURE_DIR = Path(__file__).resolve().parent
REPO_ROOT = FIGURE_DIR.parents[2]
PERFORMANCE_DIR = REPO_ROOT / "benchmark" / "performance"

GUIDE_SCALING_SOURCE = (
    PERFORMANCE_DIR / "chr22_500_guides" / "scaling_summary.tsv"
)
REFERENCE_SCALING_SOURCE = (
    PERFORMANCE_DIR
    / "chm13v2_whole_genome_20_guides"
    / "scaling_summary.tsv"
)

SOURCE_DATA_PATH = FIGURE_DIR / "figure2_source_data.tsv"
PDF_PATH = FIGURE_DIR / "figure2_performance_scaling.pdf"
PNG_PATH = FIGURE_DIR / "figure2_performance_scaling.png"
PANEL_DIR = FIGURE_DIR / "panels"

BASELINE_COLOR = "#4D4D4D"
COLUMBA_COLOR = "#0072B2"
MODE_STYLE = {
    "baseline": {"label": "Baseline", "color": BASELINE_COLOR, "marker": "o"},
    "columba": {"label": "Columba", "color": COLUMBA_COLOR, "marker": "s"},
}

SOURCE_COLUMNS = [
    "experiment",
    "panels",
    "reference",
    "reference_length_bp",
    "guide_count",
    "timing_node",
    "mode",
    "wall_seconds",
    "seconds_per_guide",
    "peak_rss_kib",
    "peak_rss_gib",
    "observed_speedup",
    "replicate_note",
    "source_file",
    "source_row",
    "source_fields",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No data rows in {path}")
    for row_number, row in enumerate(rows, start=2):
        row["_source_row"] = str(row_number)
    return rows


def single_row(rows: list[dict[str, str]], **criteria: object) -> dict[str, str]:
    matches = [
        row
        for row in rows
        if all(str(row.get(key)) == str(value) for key, value in criteria.items())
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one row for {criteria}, found {len(matches)}")
    return matches[0]


def assert_close(label: str, left: str, right: str) -> None:
    if not math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=5e-7):
        raise ValueError(f"{label}: {left} != {right}")


def environment_hostname(path: Path) -> str:
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            key, separator, value = line.rstrip("\n").partition("\t")
            if separator and key == "hostname":
                return value
    raise ValueError(f"No hostname in {path}")


def provenance_hosts(path: Path) -> str:
    hosts = sorted({row["hostname"] for row in read_tsv(path)})
    if not hosts:
        raise ValueError(f"No timing hosts in {path}")
    return ",".join(hosts)


def validate_guide_scaling(rows: list[dict[str, str]]) -> None:
    timing_20 = {
        row["mode"]: row
        for row in read_tsv(PERFORMANCE_DIR / "chr22" / "timing_summary.tsv")
        if row["mode"] in MODE_STYLE
    }
    timing_100 = {
        row["mode"]: row
        for row in read_tsv(
            PERFORMANCE_DIR / "chr22_100_guides" / "timing_summary.tsv"
        )
        if row["mode"] in MODE_STYLE
    }
    timing_500 = read_tsv(
        PERFORMANCE_DIR / "chr22_500_guides" / "timing_summary.tsv"
    )[0]

    for guide_count, timing in ((20, timing_20), (100, timing_100)):
        row = single_row(rows, reference="chr22", guide_count=guide_count)
        for mode in MODE_STYLE:
            assert_close(
                f"chr22 {guide_count} {mode} wall",
                row[f"{mode}_wall_seconds"],
                timing[mode]["median_wall_seconds"],
            )
            assert_close(
                f"chr22 {guide_count} {mode} RSS",
                row[f"{mode}_peak_rss_kib"],
                timing[mode]["median_peak_rss_kib"],
            )

    row_500 = single_row(rows, reference="chr22", guide_count=500)
    for field in (
        "baseline_wall_seconds",
        "baseline_peak_rss_kib",
        "columba_median_wall_seconds",
        "columba_median_peak_rss_kib",
    ):
        source_field = field.replace("columba_median", "columba")
        assert_close(f"chr22 500 {field}", row_500[source_field], timing_500[field])


def validate_reference_scaling(rows: list[dict[str, str]]) -> None:
    checks = [
        (
            single_row(rows, reference="chr22", guide_count=100),
            read_tsv(
                PERFORMANCE_DIR
                / "chr22_100_guides"
                / "timing_tux05_summary.tsv"
            )[0],
            "columba_median_wall_seconds",
            "columba_median_peak_rss_kib",
        ),
        (
            single_row(rows, reference="chr2", guide_count=100),
            read_tsv(PERFORMANCE_DIR / "chr2_100_guides" / "timing_summary.tsv")[0],
            "columba_median_wall_seconds",
            "columba_median_peak_rss_kib",
        ),
        (
            single_row(rows, reference="CHM13v2_whole_genome", guide_count=20),
            read_tsv(
                PERFORMANCE_DIR
                / "chm13v2_whole_genome_20_guides"
                / "timing_summary.tsv"
            )[0],
            "columba_warm_median_wall_seconds",
            "columba_warm_median_peak_rss_kib",
        ),
    ]

    for scaling, timing, columba_wall_field, columba_rss_field in checks:
        label = scaling["reference"]
        assert_close(
            f"{label} baseline wall",
            scaling["baseline_wall_seconds"],
            timing["baseline_wall_seconds"],
        )
        assert_close(
            f"{label} baseline RSS",
            scaling["baseline_peak_rss_kib"],
            timing["baseline_peak_rss_kib"],
        )
        assert_close(
            f"{label} Columba wall",
            scaling["columba_wall_seconds"],
            timing[columba_wall_field],
        )
        assert_close(
            f"{label} Columba RSS",
            scaling["columba_peak_rss_kib"],
            timing[columba_rss_field],
        )


def replicate_note(row: dict[str, str], mode: str, matched: bool) -> str:
    count = int(row[f"{mode}_replicates"])
    if matched and row["reference"] == "CHM13v2_whole_genome":
        if mode == "baseline":
            return "one controlled measured baseline replicate on tux05"
        return "median of two sequential warm-cache Columba replicates on tux05"
    method = "one measured replicate" if count == 1 else f"median of {count} measured replicates"
    return f"{method}{' on tux05' if matched else ''}"


def make_source_record(
    *,
    experiment: str,
    panels: str,
    row: dict[str, str],
    mode: str,
    timing_node: str,
    source_path: Path,
    matched: bool,
) -> dict[str, str]:
    peak_rss_kib = int(row[f"{mode}_peak_rss_kib"])
    return {
        "experiment": experiment,
        "panels": panels,
        "reference": row["reference"],
        "reference_length_bp": row["reference_length_bp"],
        "guide_count": row["guide_count"],
        "timing_node": timing_node,
        "mode": mode,
        "wall_seconds": row[f"{mode}_wall_seconds"],
        "seconds_per_guide": row[f"{mode}_seconds_per_guide"],
        "peak_rss_kib": str(peak_rss_kib),
        "peak_rss_gib": f"{peak_rss_kib / (1024 * 1024):.6f}",
        "observed_speedup": row["observed_speedup"],
        "replicate_note": replicate_note(row, mode, matched),
        "source_file": source_path.relative_to(REPO_ROOT).as_posix(),
        "source_row": row["_source_row"],
        "source_fields": ";".join(
            (
                f"{mode}_wall_seconds",
                f"{mode}_seconds_per_guide",
                f"{mode}_peak_rss_kib",
                "observed_speedup",
                f"{mode}_replicates",
            )
        ),
    }


def build_source_data() -> list[dict[str, str]]:
    guide_rows = read_tsv(GUIDE_SCALING_SOURCE)
    reference_rows = read_tsv(REFERENCE_SCALING_SOURCE)
    validate_guide_scaling(guide_rows)
    validate_reference_scaling(reference_rows)

    node_by_guide_count = {
        20: environment_hostname(PERFORMANCE_DIR / "chr22" / "environment.txt"),
        100: environment_hostname(
            PERFORMANCE_DIR / "chr22_100_guides" / "environment.txt"
        ),
        500: provenance_hosts(
            PERFORMANCE_DIR / "chr22_500_guides" / "timing_provenance.tsv"
        ),
    }

    records: list[dict[str, str]] = []
    for row in sorted(guide_rows, key=lambda item: int(item["guide_count"])):
        for mode in MODE_STYLE:
            records.append(
                make_source_record(
                    experiment="chr22_guide_count_scaling",
                    panels="A",
                    row=row,
                    mode=mode,
                    timing_node=node_by_guide_count[int(row["guide_count"])],
                    source_path=GUIDE_SCALING_SOURCE,
                    matched=False,
                )
            )

    reference_order = {"chr22": 0, "chr2": 1, "CHM13v2_whole_genome": 2}
    for row in sorted(reference_rows, key=lambda item: reference_order[item["reference"]]):
        for mode in MODE_STYLE:
            records.append(
                make_source_record(
                    experiment="matched_reference_size_scaling",
                    panels="B,C",
                    row=row,
                    mode=mode,
                    timing_node=row["timing_node"],
                    source_path=REFERENCE_SCALING_SOURCE,
                    matched=True,
                )
            )
    return records


def write_source_data(records: list[dict[str, str]]) -> None:
    with SOURCE_DATA_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def compact_log_tick(value: float, _position: float) -> str:
    if value >= 1_000:
        return f"{value / 1_000:g}k"
    if value >= 1:
        return f"{value:g}"
    return f"{value:.2g}"


def reference_tick_label(row: dict[str, str]) -> str:
    reference = (
        "CHM13v2"
        if row["reference"] == "CHM13v2_whole_genome"
        else row["reference"]
    )
    length = int(row["reference_length_bp"])
    formatted_length = (
        f"{length / 1_000_000_000:.2f} Gb"
        if length >= 1_000_000_000
        else f"{length / 1_000_000:.1f} Mb"
    )
    return f"{reference}\n{formatted_length}"


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
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def mode_legend_handles(standalone: bool = False) -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            color=MODE_STYLE[mode]["color"],
            marker=MODE_STYLE[mode]["marker"],
            markersize=6 if standalone else 5,
            markeredgecolor="white",
            markeredgewidth=0.7,
            linewidth=1.45,
            label=MODE_STYLE[mode]["label"],
        )
        for mode in MODE_STYLE
    ]


def configure_axis(axis: plt.Axes, standalone: bool = False) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#555555")
    axis.spines["bottom"].set_color("#555555")
    axis.tick_params(
        axis="both",
        which="major",
        labelsize=8.8 if standalone else 7,
        length=3,
        width=0.7,
    )
    axis.tick_params(axis="both", which="minor", length=0)
    axis.grid(axis="y", which="major", color="#E2E2E2", linewidth=0.55)
    axis.set_axisbelow(True)


def records_for(
    records: list[dict[str, str]], experiment: str, mode: str
) -> list[dict[str, str]]:
    return sorted(
        [
            record
            for record in records
            if record["experiment"] == experiment and record["mode"] == mode
        ],
        key=lambda record: int(record["reference_length_bp"])
        if experiment == "matched_reference_size_scaling"
        else int(record["guide_count"]),
    )


def plot_series(
    axis: plt.Axes,
    rows: list[dict[str, str]],
    x_field: str,
    y_field: str,
    mode: str,
) -> None:
    style = MODE_STYLE[mode]
    axis.plot(
        [float(row[x_field]) for row in rows],
        [float(row[y_field]) for row in rows],
        color=style["color"],
        linewidth=1.45,
        marker=style["marker"],
        markersize=4.8,
        markeredgecolor="white",
        markeredgewidth=0.7,
        label=style["label"],
        zorder=3,
    )


def annotate_ratios(
    axis: plt.Axes,
    baseline: list[dict[str, str]],
    columba: list[dict[str, str]],
    x_field: str,
    y_field: str,
    ratio_field: str | None = None,
    standalone: bool = False,
) -> None:
    for baseline_row, columba_row in zip(baseline, columba, strict=True):
        baseline_value = float(baseline_row[y_field])
        columba_value = float(columba_row[y_field])
        ratio = (
            float(baseline_row[ratio_field])
            if ratio_field is not None
            else columba_value / baseline_value
        )
        axis.text(
            float(baseline_row[x_field]),
            math.sqrt(baseline_value * columba_value),
            f"{ratio:.1f}x",
            ha="center",
            va="center",
            fontsize=8.2 if standalone else 6.2,
            color="#333333",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.7, "alpha": 0.9},
            zorder=4,
        )


def render_guide_count_panel(
    axis: plt.Axes, records: list[dict[str, str]], standalone: bool = False
) -> None:
    baseline = records_for(records, "chr22_guide_count_scaling", "baseline")
    columba = records_for(records, "chr22_guide_count_scaling", "columba")
    plot_series(axis, baseline, "guide_count", "wall_seconds", "baseline")
    plot_series(axis, columba, "guide_count", "wall_seconds", "columba")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(15, 650)
    axis.set_ylim(5, 30_000)
    axis.set_xticks([20, 100, 500])
    axis.set_xticklabels(["20", "100", "500"])
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.yaxis.set_major_formatter(FuncFormatter(compact_log_tick))
    axis.set_xlabel("Number of guides")
    axis.set_ylabel("Wall time (s, log scale)")
    axis.set_title("Guide-count scaling (chr22)", loc="left", pad=6)
    annotate_ratios(
        axis,
        baseline,
        columba,
        "guide_count",
        "wall_seconds",
        "observed_speedup",
        standalone,
    )
    if standalone:
        axis.legend(handles=mode_legend_handles(standalone=True), loc="upper left", frameon=False)


def render_reference_size_panel(
    axis: plt.Axes, records: list[dict[str, str]], standalone: bool = False
) -> None:
    baseline = records_for(records, "matched_reference_size_scaling", "baseline")
    columba = records_for(records, "matched_reference_size_scaling", "columba")
    labels = [reference_tick_label(row) for row in baseline]
    lengths = [float(row["reference_length_bp"]) for row in baseline]
    plot_series(
        axis, baseline, "reference_length_bp", "seconds_per_guide", "baseline"
    )
    plot_series(
        axis, columba, "reference_length_bp", "seconds_per_guide", "columba"
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(35_000_000, 4_500_000_000)
    axis.set_ylim(0.5, 3_000)
    axis.set_xticks(lengths)
    axis.set_xticklabels(labels)
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.yaxis.set_major_formatter(FuncFormatter(compact_log_tick))
    axis.set_xlabel("Reference length")
    axis.set_ylabel("Wall time per guide (s, log scale)")
    axis.set_title("Reference-size scaling (tux05)", loc="left", pad=6)
    annotate_ratios(
        axis,
        baseline,
        columba,
        "reference_length_bp",
        "seconds_per_guide",
        "observed_speedup",
        standalone,
    )
    if standalone:
        axis.legend(handles=mode_legend_handles(standalone=True), loc="upper left", frameon=False)


def render_memory_panel(
    axis: plt.Axes, records: list[dict[str, str]], standalone: bool = False
) -> None:
    baseline = records_for(records, "matched_reference_size_scaling", "baseline")
    columba = records_for(records, "matched_reference_size_scaling", "columba")
    labels = [reference_tick_label(row) for row in baseline]
    lengths = [float(row["reference_length_bp"]) for row in baseline]
    plot_series(axis, baseline, "reference_length_bp", "peak_rss_gib", "baseline")
    plot_series(axis, columba, "reference_length_bp", "peak_rss_gib", "columba")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(35_000_000, 4_500_000_000)
    axis.set_ylim(0.08, 22)
    axis.set_xticks(lengths)
    axis.set_xticklabels(labels)
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.yaxis.set_major_formatter(FuncFormatter(compact_log_tick))
    axis.set_xlabel("Reference length")
    axis.set_ylabel("Peak RSS (GiB, log scale)")
    axis.set_title("Peak memory (tux05)", loc="left", pad=6)
    annotate_ratios(
        axis,
        baseline,
        columba,
        "reference_length_bp",
        "peak_rss_gib",
        None,
        standalone,
    )
    if standalone:
        axis.legend(handles=mode_legend_handles(standalone=True), loc="upper left", frameon=False)


def render_figure(records: list[dict[str, str]]) -> None:
    apply_style(standalone=False)

    figure, axes = plt.subplots(1, 3, figsize=(8.3, 3.6))
    figure.subplots_adjust(left=0.073, right=0.99, bottom=0.29, top=0.80, wspace=0.43)

    axis_a, axis_b, axis_c = axes
    for axis in axes:
        configure_axis(axis)

    render_guide_count_panel(axis_a, records)
    render_reference_size_panel(axis_b, records)
    render_memory_panel(axis_c, records)

    for label, axis in zip(("A", "B", "C"), axes, strict=True):
        axis.text(
            -0.17,
            1.16,
            label,
            transform=axis.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
            ha="left",
        )

    legend_handles = mode_legend_handles()
    figure.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.54, 0.985),
        ncol=2,
        frameon=False,
        handlelength=2.0,
        columnspacing=1.5,
        fontsize=7.5,
    )

    figure.text(
        0.073,
        0.045,
        "A: 20/100-guide values are medians (n=3 per mode);\n"
        "500-guide baseline n=1 and Columba median n=3.",
        ha="left",
        va="bottom",
        fontsize=6.1,
        color="#333333",
    )
    figure.text(
        0.385,
        0.045,
        "B-C: matched tux05 runs. chr22/chr2 baseline n=1 and Columba median n=3;\n"
        "whole genome baseline n=1 and Columba median of two sequential warm-cache runs.",
        ha="left",
        va="bottom",
        fontsize=6.1,
        color="#333333",
    )

    pdf_metadata = {
        "Title": "Figure 2: CRISPRapido performance scaling with Columba",
        "Subject": "Guide-count, reference-size, and memory scaling",
        "Creator": "make_figure2.py",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(PDF_PATH, format="pdf", metadata=pdf_metadata, facecolor="white")
    figure.savefig(
        PNG_PATH,
        format="png",
        dpi=320,
        metadata={"Software": "make_figure2.py"},
        facecolor="white",
    )
    plt.close(figure)


def save_standalone(
    figure: plt.Figure, name: str, subject: str
) -> tuple[Path, Path]:
    pdf_path = PANEL_DIR / f"{name}.pdf"
    png_path = PANEL_DIR / f"{name}.png"
    metadata = {
        "Title": subject,
        "Subject": subject,
        "Creator": "make_figure2.py",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(pdf_path, format="pdf", metadata=metadata, facecolor="white")
    figure.savefig(
        png_path,
        format="png",
        dpi=320,
        metadata={"Software": "make_figure2.py"},
        facecolor="white",
    )
    plt.close(figure)
    return pdf_path, png_path


def render_standalone_panel(
    records: list[dict[str, str]],
    name: str,
    subject: str,
    renderer,
    note: str,
) -> tuple[Path, Path]:
    apply_style(standalone=True)
    figure, axis = plt.subplots(figsize=(6.4, 4.6))
    figure.subplots_adjust(left=0.14, right=0.98, bottom=0.25, top=0.88)
    configure_axis(axis, standalone=True)
    renderer(axis, records, standalone=True)
    figure.text(
        0.14,
        0.045,
        note,
        ha="left",
        va="bottom",
        fontsize=7.8,
        color="#333333",
    )
    return save_standalone(figure, name, subject)


def render_standalone_panels(records: list[dict[str, str]]) -> list[Path]:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    specs = [
        (
            "panelA_guide_count_scaling",
            "Figure 2A: Guide-count scaling on chr22",
            render_guide_count_panel,
            "20/100 guides: medians of three runs per mode.\n"
            "500 guides: one baseline run and median of three Columba runs.",
        ),
        (
            "panelB_reference_size_scaling",
            "Figure 2B: Matched-node reference-size scaling",
            render_reference_size_panel,
            "Matched tux05 runs. chr22/chr2: baseline n=1, Columba median n=3.\n"
            "Whole genome: baseline n=1, median of two sequential warm-cache Columba runs.",
        ),
        (
            "panelC_memory_scaling",
            "Figure 2C: Matched-node peak memory scaling",
            render_memory_panel,
            "Matched tux05 runs. Peak aggregate RSS; ratios are Columba/baseline.\n"
            "Whole-genome Columba value uses two sequential warm-cache replicates.",
        ),
    ]
    for spec in specs:
        outputs.extend(render_standalone_panel(records, *spec))
    return outputs


def main() -> None:
    records = build_source_data()
    write_source_data(records)
    render_figure(records)
    outputs = [SOURCE_DATA_PATH, PDF_PATH, PNG_PATH]
    outputs.extend(render_standalone_panels(records))
    for output in outputs:
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {output}")
        print(f"wrote {output.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
