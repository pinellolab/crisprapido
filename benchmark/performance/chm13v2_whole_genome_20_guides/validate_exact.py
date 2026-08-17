#!/usr/bin/env python3
import argparse
import csv
import hashlib
import os
import subprocess
from collections import defaultdict
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]


def clean_env():
    return {
        key: value
        for key, value in os.environ.items()
        if key not in {"RUSTFLAGS", "LIBRARY_PATH", "LD_LIBRARY_PATH"}
    }


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify(count):
    if count == 1:
        return "unique"
    if 2 <= count <= 10:
        return "low-copy"
    if count > 10:
        return "repetitive"
    return "no-exact-hit"


def parse_exact_sam(path):
    mapped = defaultdict(set)
    exact = defaultdict(set)
    with path.open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("@"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 4:
                continue
            location = (fields[2], int(fields[3]), bool(flag & 16), fields[5])
            mapped[fields[0]].add(location)
            nm = next((tag[5:] for tag in fields[11:] if tag.startswith("NM:i:")), None)
            if fields[5] == "20M" and nm == "0":
                exact[fields[0]].add(location)
    return mapped, exact


def run(command, stdout_path, stderr_path):
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        return subprocess.run(command, stdout=stdout, stderr=stderr, env=clean_env()).returncode


def main():
    parser = argparse.ArgumentParser(description="Validate the 20-guide panel at exact-match settings.")
    parser.add_argument("--reference", default=os.environ.get("CHM13_FASTA", REPO_ROOT / "../data/real_reference/chm13v2.fa"))
    parser.add_argument("--index-prefix", default=os.environ.get("COLUMBA_INDEX", REPO_ROOT / "../results/chm13_whole_genome_index/chm13v2"))
    parser.add_argument("--columba-bin", default=os.environ.get("COLUMBA_BIN", REPO_ROOT / "../columba/build_Vanilla/columba"))
    parser.add_argument("--crisprapido-bin", default=os.environ.get("CRISPRAPIDO_BIN", REPO_ROOT / "target/release/crisprapido"))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    summary_path = SCRIPT_DIR / "validation_summary.tsv"
    if summary_path.exists() and not args.force:
        raise SystemExit("Refusing to overwrite validation_summary.tsv without --force")
    reference = Path(args.reference).resolve()
    index = Path(args.index_prefix).resolve()
    columba = Path(args.columba_bin).resolve()
    crisprapido = Path(args.crisprapido_bin).resolve()
    guides = read_tsv(SCRIPT_DIR / "guides.tsv")
    if len(guides) != 20:
        raise SystemExit(f"Expected 20 guides, found {len(guides)}")
    subprocess.run(
        ["python3", str(SCRIPT_DIR / "preflight.py"), "--reference", str(reference), "--index-prefix", str(index), "--require-ready"],
        check=True,
    )

    raw = SCRIPT_DIR / "preparation_raw" / "exact_validation"
    raw.mkdir(parents=True, exist_ok=True)
    sam = raw / "panel.columba_k0.sam"
    columba_status = run(
        [
            str(columba), "-r", str(index), "-f", str(SCRIPT_DIR / "guides.fa"),
            "-a", "all", "-m", "edit", "-e", "0", "-t", "1", "-nU", "-R",
            "-o", str(sam), "-l", str(raw / "panel.columba_k0.log"),
        ],
        raw / "panel.columba_k0.stdout.txt",
        raw / "panel.columba_k0.stderr.txt",
    )
    if columba_status != 0:
        raise SystemExit(f"Direct Columba exact validation failed with exit {columba_status}")
    mapped, exact = parse_exact_sam(sam)

    rows = []
    for guide in guides:
        guide_id = guide["guide_id"]
        paf = raw / f"{guide_id}.crisprapido_k0.paf"
        stderr = raw / f"{guide_id}.crisprapido_k0.stderr.txt"
        status = run(
            [
                str(crisprapido), "-r", str(reference), "-g", guide["guide_sequence"],
                "-p", "GG", "-m", "0", "-b", "0", "-z", "0", "-f", "0.75",
                "-t", "1", "--columba-bin", str(columba), "--columba-index", str(index),
            ],
            paf,
            stderr,
        )
        exact_count = len(exact.get(guide_id, set()))
        copy_class = classify(exact_count)
        paf_count = sum(1 for line in paf.open() if line.strip()) if paf.exists() else 0
        expected_class = guide["target_copy_class"]
        if status != 0:
            validation, notes = "fail", f"CRISPRapido exit {status}"
        elif exact_count < 1:
            validation, notes = "fail", "no direct 20M/NM=0 hit"
        elif paf_count < 1:
            validation, notes = "fail", "automatic CRISPRapido k=0 returned no PAF"
        elif copy_class != expected_class:
            validation, notes = "fail", f"expected {expected_class}, observed {copy_class}"
        else:
            validation, notes = "pass", "ok"
        rows.append(
            {
                "guide_id": guide_id,
                "canonical_chromosome": guide["canonical_chromosome"],
                "target_copy_class": expected_class,
                "observed_copy_class": copy_class,
                "columba_mapped_records": len(mapped.get(guide_id, set())),
                "columba_exact_20m_nm0_records": exact_count,
                "crisprapido_paf_records": paf_count,
                "columba_exit_status": columba_status,
                "crisprapido_exit_status": status,
                "validation_status": validation,
                "notes": notes,
            }
        )

    fields = list(rows[0])
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts = {name: sum(row["observed_copy_class"] == name for row in rows) for name in ["unique", "low-copy", "repetitive"]}
    all_pass = all(row["validation_status"] == "pass" for row in rows)
    (SCRIPT_DIR / "preparation_summary.txt").write_text(
        "\n".join(
            [
                f"reference={reference}",
                f"reference_sha256={sha256(reference)}",
                f"index_prefix={index}",
                f"selected_guides={len(rows)}",
                f"unique_guides={counts['unique']}",
                f"low_copy_guides={counts['low-copy']}",
                f"repetitive_guides={counts['repetitive']}",
                f"all_exact_validation_passed={str(all_pass).lower()}",
            ]
        )
        + "\n"
    )
    if not all_pass:
        raise SystemExit(2)
    print("exact validation passed: 20/20")


if __name__ == "__main__":
    main()

