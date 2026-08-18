#!/usr/bin/env python3
import csv
import os
import re
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
REFERENCE = Path(os.environ.get('CHR22_REFERENCE', REPO_ROOT / '../data/real_reference/chm13v2_chr22.fa')).resolve()
INDEX_PREFIX = Path(os.environ.get('COLUMBA_INDEX', REPO_ROOT / '../results/chm13_chr22_index/chm13v2_chr22')).resolve()
COLUMBA_BIN = Path(os.environ.get('COLUMBA_BIN', REPO_ROOT / '../columba/build_Vanilla/columba')).resolve()
CRISPRAPIDO_BIN = Path(os.environ.get('CRISPRAPIDO_BIN', REPO_ROOT / 'target/release/crisprapido')).resolve()
OUT_DIR = SCRIPT_DIR
RAW_DIR = OUT_DIR / 'preparation_raw'
GUIDES_DIR = OUT_DIR / 'guides'
GUIDE_LEN = 20
PAM = 'GG'
INTERVALS = 100
CHROM = '22'
CANONICAL = set('ACGT')
HOMOPOLYMER_RE = re.compile(r'(A{6,}|C{6,}|G{6,}|T{6,})')


def display_path(path):
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def read_single_fasta(path):
    name = None
    seq_parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    raise ValueError('Expected one FASTA record')
                name = line[1:].split()[0]
            else:
                seq_parts.append(line.upper())
    if name is None:
        raise ValueError('No FASTA header found')
    return name, ''.join(seq_parts)


def valid_guide(seq):
    return len(seq) == GUIDE_LEN and set(seq) <= CANONICAL and HOMOPOLYMER_RE.search(seq) is None


def parse_sam(path):
    mapped = []
    exact_nm0_20m = []
    cigars = []
    nms = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line or line.startswith('@'):
                continue
            fields = line.split('\t')
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 4:
                continue
            nm = None
            for tag in fields[11:]:
                if tag.startswith('NM:i:'):
                    nm = int(tag[5:])
                    break
            mapped.append(fields)
            cigars.append(fields[5])
            nms.append('NA' if nm is None else str(nm))
            if nm == 0 and fields[5] == '20M':
                exact_nm0_20m.append(fields)
    return mapped, exact_nm0_20m, sorted(set(cigars)), sorted(set(nms))


def run_command(cmd, stdout_path, stderr_path):
    env = {k: v for k, v in os.environ.items() if k not in {'RUSTFLAGS', 'LIBRARY_PATH', 'LD_LIBRARY_PATH'}}
    with stdout_path.open('w') as out, stderr_path.open('w') as err:
        return subprocess.run(cmd, stdout=out, stderr=err, env=env).returncode


def write_fasta(path, name, seq):
    path.write_text(f'>{name}\n{seq}\n')


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    GUIDES_DIR.mkdir(exist_ok=True)
    RAW_DIR.mkdir(exist_ok=True)

    chrom, sequence = read_single_fasta(REFERENCE)
    if chrom != CHROM:
        raise SystemExit(f'Expected FASTA header {CHROM}, found {chrom}')
    n = len(sequence)

    candidates = []
    seen = set()
    all_valid_positions = 0
    for pos in range(0, n - GUIDE_LEN - len(PAM) + 1):
        guide = sequence[pos:pos + GUIDE_LEN]
        pam = sequence[pos + GUIDE_LEN:pos + GUIDE_LEN + len(PAM)]
        if pam != PAM or not valid_guide(guide):
            continue
        all_valid_positions += 1
        if guide in seen:
            continue
        seen.add(guide)
        candidates.append({'guide_sequence': guide, 'chromosome': chrom, 'protospacer_start': pos, 'strand': '+', 'PAM': pam})

    selected = []
    for interval in range(INTERVALS):
        start = (n * interval) // INTERVALS
        end = (n * (interval + 1)) // INTERVALS
        center = (start + end) // 2
        interval_candidates = [c for c in candidates if start <= c['protospacer_start'] < end]
        if not interval_candidates:
            raise SystemExit(f'No valid candidates in interval {interval + 1}: {start}-{end}')
        best = min(interval_candidates, key=lambda c: (abs(c['protospacer_start'] - center), c['protospacer_start']))
        best = dict(best)
        best['interval_number'] = interval + 1
        best['guide_id'] = f'chr22_100_guide_{interval + 1:03d}'
        selected.append(best)

    with (OUT_DIR / 'chr22_guides.tsv').open('w', newline='') as handle:
        fieldnames = ['guide_id', 'guide_sequence', 'chromosome', 'zero_based_protospacer_start', 'one_based_protospacer_start', 'strand', 'PAM', 'interval_number']
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for item in selected:
            writer.writerow({
                'guide_id': item['guide_id'],
                'guide_sequence': item['guide_sequence'],
                'chromosome': item['chromosome'],
                'zero_based_protospacer_start': item['protospacer_start'],
                'one_based_protospacer_start': item['protospacer_start'] + 1,
                'strand': item['strand'],
                'PAM': item['PAM'],
                'interval_number': item['interval_number'],
            })

    with (OUT_DIR / 'chr22_guides.fa').open('w') as handle:
        for item in selected:
            handle.write(f'>{item["guide_id"]}\n{item["guide_sequence"]}\n')
            write_fasta(GUIDES_DIR / f'{item["guide_id"]}.fa', item['guide_id'], item['guide_sequence'])

    validation_rows = []
    for item in selected:
        guide_id = item['guide_id']
        guide_seq = item['guide_sequence']
        guide_fa = GUIDES_DIR / f'{guide_id}.fa'
        columba_sam = RAW_DIR / f'{guide_id}.columba_k0.sam'
        columba_log = RAW_DIR / f'{guide_id}.columba_k0.log'
        columba_stdout = RAW_DIR / f'{guide_id}.columba_k0.stdout.txt'
        columba_stderr = RAW_DIR / f'{guide_id}.columba_k0.stderr.txt'
        columba_cmd = [str(COLUMBA_BIN), '-r', str(INDEX_PREFIX), '-f', str(guide_fa), '-a', 'all', '-m', 'edit', '-e', '0', '-t', '1', '-o', str(columba_sam), '-l', str(columba_log)]
        columba_status = run_command(columba_cmd, columba_stdout, columba_stderr)
        mapped, exact_nm0_20m, cigars, nms = parse_sam(columba_sam) if columba_sam.exists() else ([], [], [], [])

        crispr_paf = RAW_DIR / f'{guide_id}.crisprapido_k0.paf'
        crispr_stderr = RAW_DIR / f'{guide_id}.crisprapido_k0.stderr.txt'
        crispr_cmd = [str(CRISPRAPIDO_BIN), '-r', str(REFERENCE), '-g', guide_seq, '-p', PAM, '-m', '0', '-b', '0', '-z', '0', '-f', '0.75', '-t', '1', '--columba-bin', str(COLUMBA_BIN), '--columba-index', str(INDEX_PREFIX)]
        crispr_status = run_command(crispr_cmd, crispr_paf, crispr_stderr)
        paf_count = sum(1 for _ in crispr_paf.open()) if crispr_paf.exists() else 0

        if columba_status != 0:
            status = 'fail'; notes = f'Columba exit {columba_status}'
        elif crispr_status != 0:
            status = 'fail'; notes = f'CRISPRapido exit {crispr_status}'
        elif len(exact_nm0_20m) < 1:
            status = 'fail'; notes = 'No exact Columba 20M/NM=0 mapped record'
        elif paf_count < 1:
            status = 'fail'; notes = 'Columba exact hit but CRISPRapido reported zero PAF records'
        else:
            status = 'pass'; notes = 'ok'

        hit_count = len(exact_nm0_20m)
        if hit_count == 1:
            copy_class = 'unique'
        elif 2 <= hit_count <= 10:
            copy_class = 'low-copy'
        elif hit_count > 10:
            copy_class = 'repetitive'
        else:
            copy_class = 'no-exact-hit'

        validation_rows.append({
            'guide_id': guide_id,
            'columba_mapped_record_count': len(mapped),
            'columba_exact_20m_nm0_count': len(exact_nm0_20m),
            'crisprapido_paf_count': paf_count,
            'copy_class': copy_class,
            'validation_status': status,
            'notes': notes,
            'columba_exit_status': columba_status,
            'crisprapido_exit_status': crispr_status,
            'cigars': ','.join(cigars) if cigars else 'NA',
            'nms': ','.join(nms) if nms else 'NA',
        })

    with (OUT_DIR / 'validation_summary.tsv').open('w', newline='') as handle:
        fieldnames = ['guide_id', 'columba_mapped_record_count', 'columba_exact_20m_nm0_count', 'crisprapido_paf_count', 'copy_class', 'validation_status', 'notes', 'columba_exit_status', 'crisprapido_exit_status', 'cigars', 'nms']
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(validation_rows)

    counts = {key: sum(1 for row in validation_rows if row['copy_class'] == key) for key in ['unique', 'low-copy', 'repetitive', 'no-exact-hit']}
    summary_text = f'''candidate_positions={all_valid_positions}
unique_candidate_guides={len(candidates)}
selected_guides={len(selected)}
unique_guides={counts['unique']}
low_copy_guides={counts['low-copy']}
repetitive_guides={counts['repetitive']}
no_exact_hit_guides={counts['no-exact-hit']}
all_columba_validated={all(row['columba_exact_20m_nm0_count'] >= 1 for row in validation_rows)}
all_crisprapido_validated={all(row['crisprapido_paf_count'] >= 1 and row['validation_status'] == 'pass' for row in validation_rows)}
'''
    (OUT_DIR / 'preparation_summary.txt').write_text(summary_text)
    (OUT_DIR / 'README.md').write_text(f'''# CHM13 chr22 100-Guide Performance Benchmark

This package extends the 20-guide chr22 pilot to a deterministic 100-guide panel. It compares original CRISPRapido sliding-window candidate generation with automatic Columba candidate generation followed by the same WFA2 verification, filtering, PAM extraction, CFD scoring, and PAF reporting.

## Guide Panel

Reference: `../data/real_reference/chm13v2_chr22.fa`  
Columba index prefix: `../results/chm13_chr22_index/chm13v2_chr22`  
Chromosome/header: `{chrom}`  
Reference length: `{n}` bp

Selection rules:

- Forward-strand 20-nt protospacers only.
- The two reference bases immediately following the protospacer must be `GG`.
- Guide and PAM must contain only `A/C/G/T`.
- Guides with homopolymer runs longer than 5 bases are excluded.
- Exact duplicate guide sequences are excluded before interval selection.
- The chromosome is divided into 100 approximately equal intervals.
- One candidate nearest the interval center is selected per interval.
- Ties are broken by lowest zero-based coordinate.

Preparation summary:

```text
{summary_text}```

Validation requires every guide to have at least one direct Columba `20M`/`NM:i:0` hit at edit distance 0 and at least one CRISPRapido automatic Columba PAF record at exact-match settings.

## Benchmark Parameters

- PAM: `GG`
- maximum mismatches: `0`
- maximum bulges: `1`
- maximum bulge size: `2`
- Columba candidate edit-distance bound: `candidate_e = m + b*z = 2`
- minimum match fraction: `0.75`
- threads: `1`

Correctness uses the same conservative biological-locus equivalence method as the 20-guide pilot, not raw PAF byte identity.

Raw outputs are generated under `raw/<RUN_ID>/` and excluded from Git. Preparation validation raw outputs are generated under `preparation_raw/` and excluded from Git. Compact summaries and checksum manifests are retained.
''')
    print(summary_text, end='')


if __name__ == '__main__':
    main()
