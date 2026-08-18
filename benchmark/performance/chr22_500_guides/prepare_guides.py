#!/usr/bin/env python3
import argparse
import csv
import hashlib
import os
import re
import subprocess
from pathlib import Path

CONFIGS = {
    "chr22_500_guides": {
        "reference_env": "CHR22_REFERENCE",
        "reference": "../data/real_reference/chm13v2_chr22.fa",
        "index": "../results/chm13_chr22_index/chm13v2_chr22",
        "chrom": "22",
        "intervals": 500,
        "prefix": "chr22_500_guide",
        "batch_size": 25,
    },
    "chr2_100_guides": {
        "reference_env": "CHR2_REFERENCE",
        "reference": "../data/real_reference/chm13v2_chr2.fa",
        "index": "../results/chm13_chr2_index/chm13v2_chr2",
        "chrom": "2",
        "intervals": 100,
        "prefix": "chr2_100_guide",
        "batch_size": 10,
    },
}

SCRIPT_DIR = Path(os.environ.get("BENCH_PACKAGE_DIR", Path(__file__).resolve().parent)).resolve()
REPO_ROOT = SCRIPT_DIR.parents[2]
CONFIG = CONFIGS[SCRIPT_DIR.name]
COLUMBA_BIN = Path(os.environ.get("COLUMBA_BIN", REPO_ROOT / "../columba/build_Vanilla/columba")).resolve()
CRISPRAPIDO_BIN = Path(os.environ.get("CRISPRAPIDO_BIN", REPO_ROOT / "target/release/crisprapido")).resolve()
REFERENCE = Path(os.environ.get(CONFIG["reference_env"], REPO_ROOT / CONFIG["reference"])).resolve()
INDEX_PREFIX = Path(os.environ.get("COLUMBA_INDEX", REPO_ROOT / CONFIG["index"])).resolve()
GUIDE_LEN = 20
PAM = "GG"
CANONICAL = set("ACGT")
HOMOPOLYMER_RE = re.compile(r"(A{6,}|C{6,}|G{6,}|T{6,})")



def display_path(path):
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(Path(path).resolve())

def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_single_fasta(path):
    name = None
    seq_parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise ValueError("Expected one FASTA record")
                name = line[1:].split()[0]
            else:
                seq_parts.append(line.upper())
    if name is None:
        raise ValueError("No FASTA header found")
    return name, "".join(seq_parts)


def valid_guide(seq):
    return len(seq) == GUIDE_LEN and set(seq) <= CANONICAL and HOMOPOLYMER_RE.search(seq) is None


def write_fasta(path, name, seq):
    path.write_text(f">{name}\n{seq}\n")


def run_command(cmd, stdout_path, stderr_path):
    env = {k: v for k, v in os.environ.items() if k not in {"RUSTFLAGS", "LIBRARY_PATH", "LD_LIBRARY_PATH"}}
    with stdout_path.open("w") as out, stderr_path.open("w") as err:
        return subprocess.run(cmd, stdout=out, stderr=err, env=env).returncode


def parse_sam(path):
    mapped = []
    exact = []
    if not path.exists():
        return mapped, exact
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 4:
                continue
            nm = None
            for tag in fields[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
                    break
            mapped.append(fields)
            if nm == 0 and fields[5] == "20M":
                exact.append(fields)
    return mapped, exact


def select_guides(sequence, chrom):
    candidates = []
    seen = set()
    all_valid_positions = 0
    for pos in range(0, len(sequence) - GUIDE_LEN - len(PAM) + 1):
        guide = sequence[pos : pos + GUIDE_LEN]
        pam = sequence[pos + GUIDE_LEN : pos + GUIDE_LEN + len(PAM)]
        if pam != PAM or not valid_guide(guide):
            continue
        all_valid_positions += 1
        if guide in seen:
            continue
        seen.add(guide)
        candidates.append({"guide_sequence": guide, "chromosome": chrom, "protospacer_start": pos, "strand": "+", "PAM": pam})

    selected = []
    intervals = CONFIG["intervals"]
    for interval in range(intervals):
        start = (len(sequence) * interval) // intervals
        end = (len(sequence) * (interval + 1)) // intervals
        center = (start + end) // 2
        used_guides = {item["guide_sequence"] for item in selected}
        unused = [c for c in candidates if c["guide_sequence"] not in used_guides]
        interval_candidates = [c for c in unused if start <= c["protospacer_start"] < end]
        pool = interval_candidates if interval_candidates else unused
        if not pool:
            raise SystemExit(f"No unused valid candidates remain for interval {interval + 1}: {start}-{end}")
        best = min(pool, key=lambda c: (abs(c["protospacer_start"] - center), c["protospacer_start"]))
        best = dict(best)
        best["interval_number"] = interval + 1
        best["inside_interval"] = start <= best["protospacer_start"] < end
        best["guide_id"] = f'{CONFIG["prefix"]}_{interval + 1:03d}'
        selected.append(best)
    return all_valid_positions, candidates, selected


def write_panel(selected):
    guides_dir = SCRIPT_DIR / "guides"
    guides_dir.mkdir(exist_ok=True)
    with (SCRIPT_DIR / "guides.tsv").open("w", newline="") as handle:
        fields = [
            "guide_id",
            "guide_sequence",
            "chromosome",
            "zero_based_protospacer_start",
            "one_based_protospacer_start",
            "strand",
            "PAM",
            "interval_number",
            "inside_interval",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for item in selected:
            writer.writerow(
                {
                    "guide_id": item["guide_id"],
                    "guide_sequence": item["guide_sequence"],
                    "chromosome": item["chromosome"],
                    "zero_based_protospacer_start": item["protospacer_start"],
                    "one_based_protospacer_start": item["protospacer_start"] + 1,
                    "strand": item["strand"],
                    "PAM": item["PAM"],
                    "interval_number": item["interval_number"],
                    "inside_interval": item.get("inside_interval", True),
                }
            )
    with (SCRIPT_DIR / "guides.fa").open("w") as handle:
        for item in selected:
            handle.write(f'>{item["guide_id"]}\n{item["guide_sequence"]}\n')
            write_fasta(guides_dir / f'{item["guide_id"]}.fa', item["guide_id"], item["guide_sequence"])


def validate_guides(selected):
    raw_dir = SCRIPT_DIR / "preparation_raw"
    raw_dir.mkdir(exist_ok=True)
    rows = []
    for item in selected:
        guide_id = item["guide_id"]
        guide_seq = item["guide_sequence"]
        guide_fa = SCRIPT_DIR / "guides" / f"{guide_id}.fa"
        columba_sam = raw_dir / f"{guide_id}.columba_k0.sam"
        columba_log = raw_dir / f"{guide_id}.columba_k0.log"
        columba_stdout = raw_dir / f"{guide_id}.columba_k0.stdout.txt"
        columba_stderr = raw_dir / f"{guide_id}.columba_k0.stderr.txt"
        columba_cmd = [
            str(COLUMBA_BIN),
            "-r",
            str(INDEX_PREFIX),
            "-f",
            str(guide_fa),
            "-a",
            "all",
            "-m",
            "edit",
            "-e",
            "0",
            "-t",
            "1",
            "-o",
            str(columba_sam),
            "-l",
            str(columba_log),
        ]
        columba_status = run_command(columba_cmd, columba_stdout, columba_stderr)
        mapped, exact = parse_sam(columba_sam)

        crispr_paf = raw_dir / f"{guide_id}.crisprapido_k0.paf"
        crispr_stderr = raw_dir / f"{guide_id}.crisprapido_k0.stderr.txt"
        crispr_stdout = raw_dir / f"{guide_id}.crisprapido_k0.stdout.txt"
        crispr_cmd = [
            str(CRISPRAPIDO_BIN),
            "-r",
            str(REFERENCE),
            "-g",
            guide_seq,
            "-p",
            PAM,
            "-m",
            "0",
            "-b",
            "0",
            "-z",
            "0",
            "-f",
            "0.75",
            "-t",
            "1",
            "--columba-bin",
            str(COLUMBA_BIN),
            "--columba-index",
            str(INDEX_PREFIX),
        ]
        crispr_status = run_command(crispr_cmd, crispr_paf, crispr_stderr)
        crispr_stdout.write_text("")
        paf_count = sum(1 for _ in crispr_paf.open()) if crispr_paf.exists() else 0
        hit_count = len(exact)
        if hit_count == 1:
            copy_class = "unique"
        elif 2 <= hit_count <= 10:
            copy_class = "low-copy"
        elif hit_count > 10:
            copy_class = "repetitive"
        else:
            copy_class = "no-exact-hit"
        if columba_status != 0:
            status, notes = "fail", f"Columba exit {columba_status}"
        elif crispr_status != 0:
            status, notes = "fail", f"CRISPRapido exit {crispr_status}"
        elif hit_count < 1:
            status, notes = "fail", "No direct Columba 20M/NM=0 hit"
        elif paf_count < 1:
            status, notes = "fail", "CRISPRapido automatic k=0 reported zero PAF records"
        else:
            status, notes = "pass", "ok"
        rows.append(
            {
                "guide_id": guide_id,
                "columba_mapped_record_count": len(mapped),
                "columba_exact_20m_nm0_count": hit_count,
                "crisprapido_paf_count": paf_count,
                "copy_class": copy_class,
                "validation_status": status,
                "notes": notes,
                "columba_exit_status": columba_status,
                "crisprapido_exit_status": crispr_status,
            }
        )
    with (SCRIPT_DIR / "validation_summary.tsv").open("w", newline="") as handle:
        fields = [
            "guide_id",
            "columba_mapped_record_count",
            "columba_exact_20m_nm0_count",
            "crisprapido_paf_count",
            "copy_class",
            "validation_status",
            "notes",
            "columba_exit_status",
            "crisprapido_exit_status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_batches(selected):
    batch_size = CONFIG["batch_size"]
    with (SCRIPT_DIR / "batches.tsv").open("w", newline="") as handle:
        fields = ["batch_id", "start_ordinal", "end_ordinal", "guide_count", "guide_ids", "sha256"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for i in range(0, len(selected), batch_size):
            batch = selected[i : i + batch_size]
            guide_ids = ",".join(item["guide_id"] for item in batch)
            digest = hashlib.sha256(guide_ids.encode("ascii")).hexdigest()
            writer.writerow(
                {
                    "batch_id": f"batch_{i // batch_size + 1:03d}",
                    "start_ordinal": i + 1,
                    "end_ordinal": i + len(batch),
                    "guide_count": len(batch),
                    "guide_ids": guide_ids,
                    "sha256": digest,
                }
            )


def write_summary(all_valid_positions, candidates, selected, validation_rows, chrom, sequence):
    counts = {key: sum(1 for row in validation_rows if row["copy_class"] == key) for key in ["unique", "low-copy", "repetitive", "no-exact-hit"]}
    text = "\n".join(
        [
            f"reference={display_path(REFERENCE)}",
            f"reference_sha256={sha256_file(REFERENCE)}",
            f"chromosome={chrom}",
            f"reference_length={len(sequence)}",
            f"candidate_positions={all_valid_positions}",
            f"unique_candidate_guides={len(candidates)}",
            f"selected_guides={len(selected)}",
            f"unique_guides={counts['unique']}",
            f"low_copy_guides={counts['low-copy']}",
            f"repetitive_guides={counts['repetitive']}",
            f"no_exact_hit_guides={counts['no-exact-hit']}",
            f"all_columba_validated={all(row['columba_exact_20m_nm0_count'] >= 1 for row in validation_rows)}",
            f"all_crisprapido_validated={all(row['crisprapido_paf_count'] >= 1 and row['validation_status'] == 'pass' for row in validation_rows)}",
        ]
    )
    (SCRIPT_DIR / "preparation_summary.txt").write_text(text + "\n")


def write_readme(chrom, sequence):
    (SCRIPT_DIR / "README.md").write_text(
        f"""# {SCRIPT_DIR.name} CRISPRapido + Columba Benchmark

This package prepares a deterministic guide panel and Slurm-compatible batched
benchmark for CRISPRapido baseline sliding-window candidate generation versus
automatic Columba candidate generation with WFA2 verification and CFD reporting.

Reference: `{display_path(REFERENCE)}`  
Columba index prefix: `{display_path(INDEX_PREFIX)}`  
Chromosome/header: `{chrom}`  
Reference length: `{len(sequence)}` bp

Biological parameters for correctness and timing:

- PAM: `GG`
- max mismatches: `0`
- max bulges: `1`
- max bulge size: `2`
- minimum match fraction: `0.75`
- threads: `1`
- automatic Columba candidate bound: `candidate_e = m + b*z = 2`

Guide selection is deterministic: the reference is divided into
{CONFIG["intervals"]} intervals, and one unique forward-strand 20 nt guide with
an immediately adjacent `GG` PAM is selected nearest each interval center. If a nominal interval has no valid guide after filters, the nearest unused valid guide globally is selected and marked with `inside_interval=false`.
Guides with non-ACGT bases or homopolymer runs longer than five bases are
excluded before interval selection.

Raw outputs are excluded by `.gitignore`. Use `prepare_guides.py` to regenerate
the panel and exact-match validation summaries. Use the Slurm scripts for
batched correctness and measured runs, then `aggregate_batches.py` to combine
completed batches.

The Slurm scripts intentionally do not set a partition or account. Provide
site-specific values with `SBATCH_PARTITION`, `SBATCH_ACCOUNT`, and related
environment variables when submitting.
"""
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-validation", action="store_true")
    args = parser.parse_args()
    chrom, sequence = read_single_fasta(REFERENCE)
    if chrom != CONFIG["chrom"]:
        raise SystemExit(f"Expected FASTA header {CONFIG['chrom']}, found {chrom}")
    all_valid_positions, candidates, selected = select_guides(sequence, chrom)
    write_panel(selected)
    write_batches(selected)
    validation_rows = []
    validation_path = SCRIPT_DIR / "validation_summary.tsv"
    if validation_path.exists() and not args.skip_validation:
        with validation_path.open() as handle:
            validation_rows = list(csv.DictReader(handle, delimiter="\t"))
            for row in validation_rows:
                row["columba_exact_20m_nm0_count"] = int(row["columba_exact_20m_nm0_count"])
                row["crisprapido_paf_count"] = int(row["crisprapido_paf_count"])
    elif not args.skip_validation:
        validation_rows = validate_guides(selected)
    else:
        validation_rows = [
            {
                "copy_class": "not_validated",
                "columba_exact_20m_nm0_count": 0,
                "crisprapido_paf_count": 0,
                "validation_status": "not_run",
            }
            for _ in selected
        ]
    write_summary(all_valid_positions, candidates, selected, validation_rows, chrom, sequence)
    write_readme(chrom, sequence)


if __name__ == "__main__":
    main()
