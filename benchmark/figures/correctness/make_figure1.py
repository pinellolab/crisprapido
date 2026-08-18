#!/usr/bin/env python3
"""Build manuscript Figure 1 from current source and finalized correctness TSVs."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.text import Text


FIGURE_DIR = Path(__file__).resolve().parent
REPO_ROOT = FIGURE_DIR.parents[2]
PERFORMANCE_DIR = REPO_ROOT / "benchmark" / "performance"

SOURCE_DATA_PATH = FIGURE_DIR / "figure1_source_data.tsv"
PDF_PATH = FIGURE_DIR / "figure1_correctness.pdf"
PNG_PATH = FIGURE_DIR / "figure1_correctness.png"
PANEL_DIR = FIGURE_DIR / "panels"

BASELINE_COLOR = "#4D4D4D"
COLUMBA_COLOR = "#0072B2"
LIGHT_GRAY = "#E6E6E6"
LIGHT_BLUE = "#DCEEF8"
SHARED_FILL = "#F0F5F4"
SHARED_EDGE = "#6E7D7A"

SOURCE_COLUMNS = [
    "record_type",
    "panels",
    "experiment",
    "display_label",
    "reference",
    "guide_count",
    "validation_metric",
    "baseline_raw_records",
    "columba_raw_records",
    "baseline_valid_loci",
    "columba_valid_loci",
    "baseline_recovered_loci",
    "baseline_missing_loci",
    "columba_only_valid_loci",
    "columba_invalid_records",
    "baseline_recall_percent",
    "controlled_configurations_total",
    "controlled_configurations_byte_identical",
    "workflow_lane",
    "workflow_step_order",
    "workflow_step",
    "source_file",
    "source_row",
    "source_fields",
    "interpretation",
]


WORKFLOW_STEPS = [
    {
        "lane": "original",
        "order": 1,
        "label": "Guide +\nreference FASTA",
        "source": "src/main.rs:45-52;src/main.rs:188-191",
        "tokens": (("src/main.rs", "reference: PathBuf"), ("src/main.rs", "guide: String")),
    },
    {
        "lane": "original",
        "order": 2,
        "label": "Overlapping\nsliding windows",
        "source": "src/main.rs:194-200",
        "tokens": (("src/main.rs", "step_by(step_size)"),),
    },
    {
        "lane": "original",
        "order": 3,
        "label": "Ends-free WFA2\nverification",
        "source": "src/main.rs:201-268;src/verification.rs:629-699",
        "tokens": (("src/verification.rs", "AlignmentSpan::EndsFree"),),
    },
    {
        "lane": "original",
        "order": 4,
        "label": "PAM + CRISPR\nfilters",
        "source": "src/verification.rs:179-242;src/verification.rs:559-589;src/reporting.rs:43-129",
        "tokens": (("src/verification.rs", "validation_passes"), ("src/reporting.rs", "pam_matches_requested")),
    },
    {
        "lane": "original",
        "order": 5,
        "label": "CFD scoring +\nPAF output",
        "source": "src/reporting.rs:134-285;src/cfd_score.rs:147-180",
        "tokens": (("src/reporting.rs", "get_cfd_score"), ("src/reporting.rs", "cg:Z")),
    },
    {
        "lane": "columba",
        "order": 1,
        "label": "Guide + FASTA +\nColumba index",
        "source": "src/main.rs:45-52;src/main.rs:98-104;src/main.rs:144-176",
        "tokens": (("src/main.rs", "columba_index: Option<PathBuf>"),),
    },
    {
        "lane": "columba",
        "order": 2,
        "label": "Columba candidates\n(e = m + b*z)\n+ exact-set union",
        "source": "src/columba.rs:120-129;src/columba.rs:131-229;src/columba.rs:265-285",
        "tokens": (("src/columba.rs", "candidate_edit_distance_bound"), ("src/columba.rs", "run_columba_candidate_generation")),
    },
    {
        "lane": "columba",
        "order": 3,
        "label": "SAM coordinates +\nreference span",
        "source": "src/columba.rs:288-335;src/columba.rs:338-435;src/verification.rs:320-380",
        "tokens": (("src/columba.rs", "parse_columba_sam_file"), ("src/columba.rs", "cigar_reference_span")),
    },
    {
        "lane": "columba",
        "order": 4,
        "label": "Anchored WFA2 +\nbounded fallback",
        "source": "src/verification.rs:381-449;src/verification.rs:472-556;src/verification.rs:591-627",
        "tokens": (("src/verification.rs", "anchored_candidate_alignment"), ("src/verification.rs", "fallback_intervals")),
    },
    {
        "lane": "columba",
        "order": 5,
        "label": "PAM + CRISPR\nfilters",
        "source": "src/verification.rs:179-242;src/verification.rs:303-312;src/verification.rs:559-589;src/reporting.rs:43-129",
        "tokens": (("src/verification.rs", "hit_matches_requested_pam"), ("src/reporting.rs", "preferred_overlap_hit")),
    },
    {
        "lane": "columba",
        "order": 6,
        "label": "CFD scoring +\nPAF output",
        "source": "src/reporting.rs:134-285;src/cfd_score.rs:147-180",
        "tokens": (("src/cfd_score.rs", "pub fn get_cfd_score"), ("src/reporting.rs", "cf:f")),
    },
]


CORRECTNESS_EXPERIMENTS = [
    {
        "experiment": "chr22_20_guides",
        "display_label": "chr22 - 20 guides",
        "reference": "chr22",
        "guide_count": 20,
        "path": PERFORMANCE_DIR / "chr22" / "correctness_summary.tsv",
    },
    {
        "experiment": "chr22_100_guides",
        "display_label": "chr22 - 100 guides",
        "reference": "chr22",
        "guide_count": 100,
        "path": PERFORMANCE_DIR / "chr22_100_guides" / "correctness_summary.tsv",
    },
    {
        "experiment": "chr22_500_guides",
        "display_label": "chr22 - 500 guides",
        "reference": "chr22",
        "guide_count": 500,
        "path": PERFORMANCE_DIR / "chr22_500_guides" / "correctness_summary.tsv",
    },
    {
        "experiment": "chr2_100_guides",
        "display_label": "chr2 - 100 guides",
        "reference": "chr2",
        "guide_count": 100,
        "path": PERFORMANCE_DIR / "chr2_100_guides" / "correctness_summary.tsv",
    },
    {
        "experiment": "chm13v2_whole_genome_20_guides",
        "display_label": "CHM13v2 - 20 guides",
        "reference": "CHM13v2_whole_genome",
        "guide_count": 20,
        "path": PERFORMANCE_DIR
        / "chm13v2_whole_genome_20_guides"
        / "correctness_summary.tsv",
    },
]


def empty_source_record() -> dict[str, str]:
    return {column: "" for column in SOURCE_COLUMNS}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No data rows in {path}")
    for row_number, row in enumerate(rows, start=2):
        row["_source_row"] = str(row_number)
    return rows


def validate_workflow_sources() -> None:
    source_cache: dict[str, str] = {}
    for step in WORKFLOW_STEPS:
        for relative_path, token in step["tokens"]:
            if relative_path not in source_cache:
                source_cache[relative_path] = (REPO_ROOT / relative_path).read_text(
                    encoding="utf-8"
                )
            if token not in source_cache[relative_path]:
                raise ValueError(
                    f"Workflow token {token!r} not found in {relative_path}; "
                    "review Figure 1 against the current implementation"
                )


def workflow_source_records() -> list[dict[str, str]]:
    validate_workflow_sources()
    records = []
    for step in WORKFLOW_STEPS:
        record = empty_source_record()
        record.update(
            {
                "record_type": "workflow_step",
                "panels": "A",
                "experiment": "integration_workflow",
                "validation_metric": "source_audited_implementation_flow",
                "workflow_lane": step["lane"],
                "workflow_step_order": str(step["order"]),
                "workflow_step": step["label"].replace("\n", " "),
                "source_file": ";".join(
                    sorted({location.split(":", 1)[0] for location in step["source"].split(";")})
                ),
                "source_row": step["source"],
                "source_fields": ";".join(token for _, token in step["tokens"]),
                "interpretation": "Current production workflow step.",
            }
        )
        records.append(record)
    return records


def correctness_source_records() -> list[dict[str, str]]:
    records = []
    for experiment in CORRECTNESS_EXPERIMENTS:
        row = read_tsv(experiment["path"])[0]
        baseline_valid = int(row["baseline_valid_loci"])
        columba_valid = int(row["columba_valid_loci"])
        recovered = int(row["shared_baseline_loci"])
        missing = int(row["baseline_missing_from_columba"])
        columba_only = int(row["columba_only_valid_loci"])
        invalid = int(row["columba_invalid_records"])

        if baseline_valid != recovered + missing:
            raise ValueError(f"Inconsistent baseline locus accounting in {experiment['path']}")
        if columba_valid != recovered + columba_only:
            raise ValueError(f"Inconsistent Columba locus accounting in {experiment['path']}")
        recall = 100.0 * recovered / baseline_valid if baseline_valid else 0.0

        record = empty_source_record()
        record.update(
            {
                "record_type": "correctness_result",
                "panels": "B,C"
                if experiment["reference"] == "CHM13v2_whole_genome"
                else "C",
                "experiment": experiment["experiment"],
                "display_label": experiment["display_label"],
                "reference": experiment["reference"],
                "guide_count": str(experiment["guide_count"]),
                "validation_metric": "oracle_valid_baseline_locus_recall",
                "baseline_raw_records": row["baseline_raw_records"],
                "columba_raw_records": row["columba_raw_records"],
                "baseline_valid_loci": str(baseline_valid),
                "columba_valid_loci": str(columba_valid),
                "baseline_recovered_loci": str(recovered),
                "baseline_missing_loci": str(missing),
                "columba_only_valid_loci": str(columba_only),
                "columba_invalid_records": str(invalid),
                "baseline_recall_percent": f"{recall:.6f}",
                "source_file": experiment["path"].relative_to(REPO_ROOT).as_posix(),
                "source_row": row["_source_row"],
                "source_fields": ";".join(
                    (
                        "baseline_raw_records",
                        "columba_raw_records",
                        "baseline_valid_loci",
                        "columba_valid_loci",
                        "shared_baseline_loci",
                        "baseline_missing_from_columba",
                        "columba_only_valid_loci",
                        "columba_invalid_records",
                    )
                ),
                "interpretation": (
                    f"{recovered}/{baseline_valid} oracle-valid baseline loci recovered; "
                    f"{missing} missing; {invalid} independently invalid Columba records."
                ),
            }
        )
        records.append(record)
    return records


def controlled_source_record() -> dict[str, str]:
    path = REPO_ROOT / "benchmark" / "controlled" / "summary.tsv"
    rows = read_tsv(path)
    identical = [
        row
        for row in rows
        if row["paf_byte_identical"] == "yes"
        and row["manual_exit"] == "0"
        and row["automatic_exit"] == "0"
    ]
    if len(identical) != len(rows):
        raise ValueError("Not all controlled configurations are successful and byte-identical")

    record = empty_source_record()
    record.update(
        {
            "record_type": "controlled_result",
            "panels": "C",
            "experiment": "controlled_synthetic",
            "display_label": "Controlled synthetic",
            "validation_metric": "candidate_matched_raw_paf_identity",
            "controlled_configurations_total": str(len(rows)),
            "controlled_configurations_byte_identical": str(len(identical)),
            "source_file": path.relative_to(REPO_ROOT).as_posix(),
            "source_row": f"2-{len(rows) + 1}",
            "source_fields": "manual_exit;automatic_exit;paf_byte_identical",
            "interpretation": (
                f"{len(identical)}/{len(rows)} candidate-matched configurations produced "
                "byte-identical manual and automatic PAF; this is not an oracle locus count."
            ),
        }
    )
    return record


def build_source_data() -> list[dict[str, str]]:
    return workflow_source_records() + correctness_source_records() + [
        controlled_source_record()
    ]


def write_source_data(records: list[dict[str, str]]) -> None:
    with SOURCE_DATA_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def draw_workflow_row(
    axis: plt.Axes,
    *,
    labels: list[str],
    y: float,
    fills: list[str],
    edges: list[str],
) -> None:
    start = 0.155
    end = 0.99
    gap = 0.014
    width = (end - start - gap * (len(labels) - 1)) / len(labels)
    height = 0.23

    for index, (label, fill, edge) in enumerate(zip(labels, fills, edges, strict=True)):
        x = start + index * (width + gap)
        node = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.004,rounding_size=0.012",
            linewidth=1.0,
            edgecolor=edge,
            facecolor=fill,
            zorder=2,
        )
        axis.add_patch(node)
        axis.text(
            x + width / 2,
            y + height / 2,
            label,
            ha="center",
            va="center",
            fontsize=5.7,
            color="#222222",
            linespacing=1.15,
            zorder=3,
        )
        if index + 1 < len(labels):
            axis.add_patch(
                FancyArrowPatch(
                    (x + width + 0.002, y + height / 2),
                    (x + width + gap - 0.002, y + height / 2),
                    arrowstyle="-|>",
                    mutation_scale=7,
                    linewidth=0.8,
                    color="#707070",
                    zorder=1,
                )
            )


def configure_quantitative_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#555555")
    axis.spines["bottom"].set_color("#555555")
    axis.tick_params(axis="both", which="major", labelsize=7, length=3, width=0.7)
    axis.grid(axis="x", color="#E2E2E2", linewidth=0.55)
    axis.set_axisbelow(True)


def scale_panel_typography(axis: plt.Axes, factor: float = 1.35) -> None:
    for item in axis.findobj(match=Text):
        item.set_fontsize(item.get_fontsize() * factor)


def render_workflow(axis: plt.Axes) -> None:
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    axis.set_title("Integration workflow", loc="left", pad=5)

    original = sorted(
        [step for step in WORKFLOW_STEPS if step["lane"] == "original"],
        key=lambda step: step["order"],
    )
    columba = sorted(
        [step for step in WORKFLOW_STEPS if step["lane"] == "columba"],
        key=lambda step: step["order"],
    )

    axis.text(
        0.006,
        0.725,
        "Original\nCRISPRapido",
        ha="left",
        va="center",
        fontsize=7.2,
        fontweight="semibold",
        color=BASELINE_COLOR,
    )
    axis.text(
        0.006,
        0.295,
        "Columba-enabled\nCRISPRapido",
        ha="left",
        va="center",
        fontsize=7.2,
        fontweight="semibold",
        color=COLUMBA_COLOR,
    )

    draw_workflow_row(
        axis,
        labels=[step["label"] for step in original],
        y=0.61,
        fills=["#FAFAFA", LIGHT_GRAY, SHARED_FILL, SHARED_FILL, SHARED_FILL],
        edges=["#888888", BASELINE_COLOR, SHARED_EDGE, SHARED_EDGE, SHARED_EDGE],
    )
    draw_workflow_row(
        axis,
        labels=[step["label"] for step in columba],
        y=0.18,
        fills=["#FAFAFA", LIGHT_BLUE, "#F2F7FA", SHARED_FILL, SHARED_FILL, SHARED_FILL],
        edges=["#888888", COLUMBA_COLOR, COLUMBA_COLOR, SHARED_EDGE, SHARED_EDGE, SHARED_EDGE],
    )

    axis.text(
        0.57,
        0.035,
        "Columba replaces candidate generation; WFA2 verification, PAM/filtering, "
        "CFD scoring, and PAF reporting remain downstream.",
        ha="center",
        va="center",
        fontsize=6.4,
        color="#333333",
    )


def render_whole_genome_correctness(axis: plt.Axes, row: dict[str, str]) -> None:
    configure_quantitative_axis(axis)
    baseline_valid = int(row["baseline_valid_loci"])
    recovered = int(row["baseline_recovered_loci"])
    additional = int(row["columba_only_valid_loci"])
    columba_valid = int(row["columba_valid_loci"])
    missing = int(row["baseline_missing_loci"])
    invalid = int(row["columba_invalid_records"])
    recall = float(row["baseline_recall_percent"])

    axis.barh(1, baseline_valid, height=0.42, color=BASELINE_COLOR, zorder=3)
    axis.barh(0, recovered, height=0.42, color=BASELINE_COLOR, zorder=3)
    axis.barh(
        0,
        additional,
        left=recovered,
        height=0.42,
        color=COLUMBA_COLOR,
        zorder=3,
    )
    axis.set_yticks([1, 0])
    axis.set_yticklabels(["Baseline valid", "Columba valid"])
    axis.set_xlim(0, columba_valid * 1.08)
    axis.set_ylim(-0.65, 1.72)
    axis.set_xlabel("Independently validated loci")
    axis.set_title("Whole-genome correctness", loc="left", pad=6)

    axis.text(
        baseline_valid / 2,
        1,
        f"{baseline_valid}",
        ha="center",
        va="center",
        color="white",
        fontsize=7,
        fontweight="semibold",
    )
    axis.text(
        recovered / 2,
        0,
        f"{recovered} shared",
        ha="center",
        va="center",
        color="white",
        fontsize=6.7,
        fontweight="semibold",
    )
    axis.text(
        recovered + additional / 2,
        0,
        f"+{additional} additional\noracle-valid candidates",
        ha="center",
        va="center",
        color="white",
        fontsize=6.5,
        fontweight="semibold",
        linespacing=1.1,
    )
    axis.text(
        0.98,
        0.97,
        f"{recall:.0f}% baseline recovery\n{recovered}/{baseline_valid}; {missing} missing",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=7.2,
        fontweight="semibold",
        color=COLUMBA_COLOR,
    )
    axis.text(
        0.01,
        0.035,
        f"{invalid} independently invalid Columba records\n"
        "Validated loci (not raw PAF records)",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.1,
        color="#444444",
    )


def render_cross_scale_correctness(
    axis: plt.Axes,
    rows: list[dict[str, str]],
    controlled: dict[str, str],
) -> None:
    configure_quantitative_axis(axis)
    y_positions = list(range(len(rows)))
    recalls = [float(row["baseline_recall_percent"]) for row in rows]
    axis.barh(y_positions, [100] * len(rows), height=0.52, color=LIGHT_GRAY, zorder=1)
    axis.barh(y_positions, recalls, height=0.52, color=COLUMBA_COLOR, zorder=2)
    axis.set_yticks(y_positions)
    axis.set_yticklabels([row["display_label"] for row in rows])
    axis.invert_yaxis()
    axis.set_xlim(0, 139)
    axis.set_xticks([0, 50, 100])
    axis.set_xlabel("Recovered baseline-valid loci (%)")
    axis.set_title("Correctness across benchmark scales", loc="left", pad=6)

    for position, row in zip(y_positions, rows, strict=True):
        recovered = int(row["baseline_recovered_loci"])
        baseline_valid = int(row["baseline_valid_loci"])
        missing = int(row["baseline_missing_loci"])
        invalid = int(row["columba_invalid_records"])
        axis.text(
            97.5,
            position,
            f"{recovered}/{baseline_valid}",
            ha="right",
            va="center",
            color="white",
            fontsize=6.5,
            fontweight="semibold",
        )
        axis.text(
            102,
            position,
            f"missing {missing} | invalid {invalid}",
            ha="left",
            va="center",
            fontsize=6.2,
            color="#333333",
        )

    total = controlled["controlled_configurations_total"]
    identical = controlled["controlled_configurations_byte_identical"]
    axis.text(
        0.0,
        -0.22,
        f"Controlled synthetic: {identical}/{total} candidate-matched configurations\n"
        "produced byte-identical PAF (separate raw-output test).",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=6.1,
        color="#444444",
    )


def render_figure(records: list[dict[str, str]]) -> None:
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

    figure = plt.figure(figsize=(8.3, 5.5))
    grid = figure.add_gridspec(
        2,
        2,
        height_ratios=(0.95, 1.2),
        left=0.105,
        right=0.99,
        bottom=0.14,
        top=0.94,
        hspace=0.47,
        wspace=0.36,
    )
    axis_a = figure.add_subplot(grid[0, :])
    axis_b = figure.add_subplot(grid[1, 0])
    axis_c = figure.add_subplot(grid[1, 1])

    correctness = [
        row for row in records if row["record_type"] == "correctness_result"
    ]
    controlled = next(
        row for row in records if row["record_type"] == "controlled_result"
    )
    whole_genome = next(
        row
        for row in correctness
        if row["reference"] == "CHM13v2_whole_genome"
    )

    render_workflow(axis_a)
    render_whole_genome_correctness(axis_b, whole_genome)
    render_cross_scale_correctness(axis_c, correctness, controlled)

    for label, axis, x_position in (
        ("A", axis_a, -0.025),
        ("B", axis_b, -0.16),
        ("C", axis_c, -0.16),
    ):
        axis.text(
            x_position,
            1.10,
            label,
            transform=axis.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
            ha="left",
        )

    pdf_metadata = {
        "Title": "Figure 1: Columba integration and correctness",
        "Subject": "CRISPRapido workflow and validated locus recovery",
        "Creator": "make_figure1.py",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(PDF_PATH, format="pdf", metadata=pdf_metadata, facecolor="white")
    figure.savefig(
        PNG_PATH,
        format="png",
        dpi=320,
        metadata={"Software": "make_figure1.py"},
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
        "Creator": "make_figure1.py",
        "CreationDate": None,
        "ModDate": None,
    }
    figure.savefig(pdf_path, format="pdf", metadata=metadata, facecolor="white")
    figure.savefig(
        png_path,
        format="png",
        dpi=320,
        metadata={"Software": "make_figure1.py"},
        facecolor="white",
    )
    plt.close(figure)
    return pdf_path, png_path


def render_standalone_panels(records: list[dict[str, str]]) -> list[Path]:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    correctness = [
        row for row in records if row["record_type"] == "correctness_result"
    ]
    controlled = next(
        row for row in records if row["record_type"] == "controlled_result"
    )
    whole_genome = next(
        row
        for row in correctness
        if row["reference"] == "CHM13v2_whole_genome"
    )
    outputs: list[Path] = []

    figure_a, axis_a = plt.subplots(figsize=(10.0, 4.2))
    figure_a.subplots_adjust(left=0.035, right=0.99, bottom=0.08, top=0.90)
    render_workflow(axis_a)
    scale_panel_typography(axis_a, 1.35)
    outputs.extend(
        save_standalone(
            figure_a,
            "panelA_workflow",
            "Figure 1A: CRISPRapido and Columba-enabled workflows",
        )
    )

    figure_b, axis_b = plt.subplots(figsize=(6.4, 4.3))
    figure_b.subplots_adjust(left=0.18, right=0.98, bottom=0.16, top=0.90)
    render_whole_genome_correctness(axis_b, whole_genome)
    scale_panel_typography(axis_b, 1.30)
    outputs.extend(
        save_standalone(
            figure_b,
            "panelB_whole_genome_recovery",
            "Figure 1B: Whole-genome validated locus recovery",
        )
    )

    figure_c, axis_c = plt.subplots(figsize=(7.2, 4.7))
    figure_c.subplots_adjust(left=0.23, right=0.98, bottom=0.24, top=0.91)
    render_cross_scale_correctness(axis_c, correctness, controlled)
    scale_panel_typography(axis_c, 1.30)
    outputs.extend(
        save_standalone(
            figure_c,
            "panelC_recovery_across_benchmarks",
            "Figure 1C: Correctness across benchmark scales",
        )
    )
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
