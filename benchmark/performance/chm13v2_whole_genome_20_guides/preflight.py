#!/usr/bin/env python3
import argparse
import hashlib
import os
import shutil
from pathlib import Path

from prepare_guides import canonical_chromosome


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
INDEX_SUFFIXES = [
    ".brt",
    ".bwt",
    ".cct",
    ".fsid",
    ".headerSN.bin",
    ".meta",
    ".pos",
    ".rev.brt",
    ".sa.4",
    ".sa.bv.4",
    ".sna",
    ".txt.bin",
]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_inventory(path):
    records = []
    name = None
    length = 0
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, length))
                name = line[1:].split()[0]
                length = 0
            else:
                if name is None:
                    raise ValueError("sequence before FASTA header")
                length += len(line)
    if name is not None:
        records.append((name, length))
    return records


def main():
    parser = argparse.ArgumentParser(description="Audit whole-genome benchmark prerequisites.")
    parser.add_argument(
        "--reference",
        default=os.environ.get("CHM13_FASTA", REPO_ROOT / "../data/real_reference/chm13v2.fa"),
    )
    parser.add_argument(
        "--index-prefix",
        default=os.environ.get(
            "COLUMBA_INDEX", REPO_ROOT / "../results/chm13_whole_genome_index/chm13v2"
        ),
    )
    parser.add_argument("--output", default=SCRIPT_DIR / "reference_inventory.tsv")
    parser.add_argument("--require-ready", action="store_true")
    args = parser.parse_args()

    reference = Path(args.reference).resolve()
    index = Path(args.index_prefix).resolve()
    output = Path(args.output)
    reference_exists = reference.is_file() and os.access(reference, os.R_OK)
    records = fasta_inventory(reference) if reference_exists else []
    missing_index = [str(index) + suffix for suffix in INDEX_SUFFIXES if not Path(str(index) + suffix).is_file()]
    present_index = [Path(str(index) + suffix) for suffix in INDEX_SUFFIXES if Path(str(index) + suffix).is_file()]
    disk = shutil.disk_usage(reference.parent if reference.parent.exists() else REPO_ROOT)

    rows = [
        ("reference_path", str(reference)),
        ("reference_exists", str(reference_exists).lower()),
        ("reference_size_bytes", reference.stat().st_size if reference_exists else 0),
        ("reference_sha256", sha256(reference) if reference_exists else "missing"),
        ("fasta_records", len(records)),
        ("reference_bases", sum(length for _, length in records)),
        ("reference_headers", ",".join(name for name, _ in records) if records else "missing"),
        (
            "normalized_reference_headers",
            ",".join(
                f"{name}={canonical_chromosome(name) or 'unrecognized'}" for name, _ in records
            )
            if records
            else "missing",
        ),
        ("index_prefix", str(index)),
        ("index_complete", str(not missing_index).lower()),
        ("index_size_bytes", sum(path.stat().st_size for path in present_index)),
        ("missing_index_files", ",".join(missing_index) if missing_index else "none"),
        ("disk_available_bytes", disk.free),
        ("ready", str(reference_exists and not missing_index).lower()),
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        handle.write("metric\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")
    for key, value in rows:
        print(f"{key}={value}")
    if args.require_ready and (not reference_exists or missing_index):
        raise SystemExit(2)


if __name__ == "__main__":
    main()

