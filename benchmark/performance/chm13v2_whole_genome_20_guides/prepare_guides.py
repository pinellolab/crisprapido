#!/usr/bin/env python3
import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
GUIDE_LEN = 20
PAM = "GG"
CANONICAL = set("ACGT")
HOMOPOLYMER = re.compile(r"(A{6,}|C{6,}|G{6,}|T{6,})")
ACCESSION_TO_CHROM = {
    **{f"CP{68278 - chrom:06d}": str(chrom) for chrom in range(1, 23)},
    "CP086569": "X",
    "CP086568": "Y",
}


def clean_env():
    return {
        key: value
        for key, value in os.environ.items()
        if key not in {"RUSTFLAGS", "LIBRARY_PATH", "LD_LIBRARY_PATH"}
    }


def canonical_chromosome(header):
    token = header.split()[0]
    contig = token.rsplit("#", 1)[-1]
    plain = contig[3:] if contig.lower().startswith("chr") else contig
    canonical = plain.upper()
    if canonical in {str(value) for value in range(1, 23)} | {"X", "Y", "M"}:
        return canonical
    accession = token.split(".")[0].upper()
    return ACCESSION_TO_CHROM.get(accession)


def fasta_records(path):
    name = None
    parts = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                if name is None:
                    raise ValueError("sequence before FASTA header")
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def valid_guide(sequence):
    return (
        len(sequence) == GUIDE_LEN
        and set(sequence) <= CANONICAL
        and HOMOPOLYMER.search(sequence) is None
    )


def candidate_positions(sequence):
    last = len(sequence) - GUIDE_LEN - len(PAM)
    center = max(0, last // 2)
    for distance in range(last + 1):
        positions = sorted({center - distance, center + distance})
        for position in positions:
            if 0 <= position <= last:
                guide = sequence[position : position + GUIDE_LEN]
                pam = sequence[position + GUIDE_LEN : position + GUIDE_LEN + len(PAM)]
                if pam == PAM and valid_guide(guide):
                    yield position, guide, distance


def read_targets(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def classify(count):
    if count == 1:
        return "unique"
    if 2 <= count <= 10:
        return "low-copy"
    if count > 10:
        return "repetitive"
    return "no-exact-hit"


def parse_exact_counts(path):
    locations = {}
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
            nm = next((tag[5:] for tag in fields[11:] if tag.startswith("NM:i:")), None)
            if fields[5] != "20M" or nm != "0":
                continue
            key = (fields[2], int(fields[3]), bool(flag & 16), fields[5])
            locations.setdefault(fields[0], set()).add(key)
    return locations


def run_columba_pool(columba, index, candidates, canonical, pool_number, raw_dir):
    pool_prefix = raw_dir / f"chr{canonical}.pool_{pool_number:03d}"
    query_path = pool_prefix.with_suffix(".fa")
    sam_path = pool_prefix.with_suffix(".sam")
    log_path = pool_prefix.with_suffix(".log")
    stdout_path = pool_prefix.with_suffix(".stdout.txt")
    stderr_path = pool_prefix.with_suffix(".stderr.txt")
    with query_path.open("w") as handle:
        for candidate in candidates:
            handle.write(f">{candidate['query_name']}\n{candidate['guide_sequence']}\n")
    command = [
        str(columba), "-r", str(index), "-f", str(query_path), "-a", "all",
        "-m", "edit", "-e", "0", "-t", "1", "-nU", "-R",
        "-o", str(sam_path), "-l", str(log_path),
    ]
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, env=clean_env())
    if result.returncode != 0:
        raise RuntimeError(f"Columba pool failed for chromosome {canonical}: exit {result.returncode}")
    return parse_exact_counts(sam_path)


def select_for_target(sequence, header, target, columba, index, selected_sequences, raw_dir, pool_size, max_pools):
    canonical = target["canonical_chromosome"]
    desired_class = target["target_copy_class"]
    positions = candidate_positions(sequence)
    seen_candidates = set()
    examined = 0
    for pool_number in range(1, max_pools + 1):
        pool = []
        while len(pool) < pool_size:
            try:
                position, guide, distance = next(positions)
            except StopIteration:
                break
            if guide in selected_sequences or guide in seen_candidates:
                continue
            seen_candidates.add(guide)
            examined += 1
            pool.append(
                {
                    "query_name": f"candidate_{canonical}_{examined:06d}",
                    "guide_sequence": guide,
                    "position": position,
                    "distance": distance,
                }
            )
        if not pool:
            break
        exact = run_columba_pool(columba, index, pool, canonical, pool_number, raw_dir)
        for candidate in pool:
            count = len(exact.get(candidate["query_name"], set()))
            if classify(count) == desired_class:
                return {
                    **candidate,
                    "chromosome": header,
                    "canonical_chromosome": canonical,
                    "target_copy_class": desired_class,
                    "exact_hit_count_at_selection": count,
                    "selection_pool": pool_number,
                    "selection_candidates_examined": examined,
                }
    raise RuntimeError(
        f"No {desired_class} guide found for chromosome {canonical} after {examined} candidates"
    )


def write_panel(selected):
    guides_dir = SCRIPT_DIR / "guides"
    guides_dir.mkdir(exist_ok=False)
    fields = [
        "guide_id", "guide_sequence", "chromosome", "canonical_chromosome",
        "zero_based_protospacer_start", "one_based_protospacer_start", "strand",
        "PAM", "target_copy_class", "exact_hit_count_at_selection",
        "selection_distance_from_center", "selection_pool", "selection_candidates_examined",
    ]
    with (SCRIPT_DIR / "guides.tsv").open("w", newline="") as handle, (SCRIPT_DIR / "guides.fa").open("w") as fasta:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for item in selected:
            row = {
                "guide_id": item["guide_id"],
                "guide_sequence": item["guide_sequence"],
                "chromosome": item["chromosome"],
                "canonical_chromosome": item["canonical_chromosome"],
                "zero_based_protospacer_start": item["position"],
                "one_based_protospacer_start": item["position"] + 1,
                "strand": "+",
                "PAM": PAM,
                "target_copy_class": item["target_copy_class"],
                "exact_hit_count_at_selection": item["exact_hit_count_at_selection"],
                "selection_distance_from_center": item["distance"],
                "selection_pool": item["selection_pool"],
                "selection_candidates_examined": item["selection_candidates_examined"],
            }
            writer.writerow(row)
            fasta.write(f">{item['guide_id']}\n{item['guide_sequence']}\n")
            (guides_dir / f"{item['guide_id']}.fa").write_text(
                f">{item['guide_id']}\n{item['guide_sequence']}\n"
            )

    with (SCRIPT_DIR / "batches.tsv").open("w", newline="") as handle:
        fields = ["batch_id", "start_ordinal", "end_ordinal", "guide_count", "guide_ids", "sha256"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for offset in range(0, len(selected), 2):
            batch = selected[offset : offset + 2]
            guide_ids = ",".join(item["guide_id"] for item in batch)
            writer.writerow(
                {
                    "batch_id": f"batch_{offset // 2 + 1:03d}",
                    "start_ordinal": offset + 1,
                    "end_ordinal": offset + len(batch),
                    "guide_count": len(batch),
                    "guide_ids": guide_ids,
                    "sha256": hashlib.sha256(guide_ids.encode("ascii")).hexdigest(),
                }
            )


def run_self_test():
    assert canonical_chromosome("chr22") == "22"
    assert canonical_chromosome("2") == "2"
    assert canonical_chromosome("CP068277.2") == "1"
    assert canonical_chromosome("CP086569.2") == "X"
    assert canonical_chromosome("CHM13#0#chr1") == "1"
    assert canonical_chromosome("CHM13#0#chr22") == "22"
    assert canonical_chromosome("CHM13#0#chrX") == "X"
    assert canonical_chromosome("CHM13#0#chrM") == "M"
    assert valid_guide("ACGT" * 5)
    assert not valid_guide("AAAAAA" + "CGTACGTACGTACG")
    sequence = "A" * 20 + "GG" + "C" * 20 + "GG"
    candidates = list(candidate_positions(sequence))
    assert all(sequence[pos + 20 : pos + 22] == "GG" for pos, _, _ in candidates)
    assert classify(1) == "unique" and classify(2) == "low-copy" and classify(11) == "repetitive"
    print("self-test ok")


def main():
    parser = argparse.ArgumentParser(description="Prepare deterministic CHM13v2 whole-genome guides.")
    parser.add_argument("--reference", default=os.environ.get("CHM13_FASTA", REPO_ROOT / "../data/real_reference/chm13v2.fa"))
    parser.add_argument("--index-prefix", default=os.environ.get("COLUMBA_INDEX", REPO_ROOT / "../results/chm13_whole_genome_index/chm13v2"))
    parser.add_argument("--columba-bin", default=os.environ.get("COLUMBA_BIN", REPO_ROOT / "../columba/build_Vanilla/columba"))
    parser.add_argument("--pool-size", type=int, default=200)
    parser.add_argument("--max-pools", type=int, default=20)
    parser.add_argument("--skip-final-validation", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return

    reference = Path(args.reference).resolve()
    index = Path(args.index_prefix).resolve()
    columba = Path(args.columba_bin).resolve()
    if (SCRIPT_DIR / "guides.tsv").exists() or (SCRIPT_DIR / "guides").exists():
        raise SystemExit("Refusing to overwrite an existing panel; move it explicitly before rerunning")
    if not reference.is_file():
        raise SystemExit(f"Missing full CHM13v2 FASTA: {reference}")
    if not columba.is_file():
        raise SystemExit(f"Missing Columba executable: {columba}")
    subprocess.run(
        [sys.executable, str(SCRIPT_DIR / "preflight.py"), "--reference", str(reference), "--index-prefix", str(index), "--require-ready"],
        check=True,
    )

    targets = read_targets(SCRIPT_DIR / "panel_targets.tsv")
    target_by_chrom = {row["canonical_chromosome"]: row for row in targets}
    found = {}
    raw_dir = SCRIPT_DIR / "preparation_raw"
    raw_dir.mkdir(exist_ok=True)
    selected_sequences = set()
    for header, sequence in fasta_records(reference):
        canonical = canonical_chromosome(header)
        if canonical not in target_by_chrom or canonical in found:
            continue
        selected = select_for_target(
            sequence, header, target_by_chrom[canonical], columba, index,
            selected_sequences, raw_dir, args.pool_size, args.max_pools,
        )
        selected_sequences.add(selected["guide_sequence"])
        found[canonical] = selected
        print(
            f"selected chr{canonical} {selected['target_copy_class']} "
            f"{selected['guide_sequence']} at {selected['position']}"
        )
    missing = [row["canonical_chromosome"] for row in targets if row["canonical_chromosome"] not in found]
    if missing:
        raise SystemExit(f"Missing target chromosomes in FASTA: {','.join(missing)}")

    selected = []
    for target in targets:
        item = found[target["canonical_chromosome"]]
        item["guide_id"] = f"chm13_wg20_guide_{int(target['panel_ordinal']):03d}"
        selected.append(item)
    write_panel(selected)
    if not args.skip_final_validation:
        subprocess.run(
            [sys.executable, str(SCRIPT_DIR / "validate_exact.py"), "--reference", str(reference), "--index-prefix", str(index), "--columba-bin", str(columba)],
            check=True,
            env=clean_env(),
        )


if __name__ == "__main__":
    main()

