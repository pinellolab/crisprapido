#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BENCH_DIR="${SCRIPT_DIR}"
RAW_ROOT="${BENCH_DIR}/raw"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_DIR="${RAW_ROOT}/${RUN_ID}"

CRISPRAPIDO_BIN="${CRISPRAPIDO_BIN:-${REPO_ROOT}/target/release/crisprapido}"
COLUMBA_BIN="${COLUMBA_BIN:-${REPO_ROOT}/../columba/build_Vanilla/columba}"
REFERENCE="${CHR22_REFERENCE:-${REPO_ROOT}/../data/real_reference/chm13v2_chr22.fa}"
COLUMBA_INDEX="${COLUMBA_INDEX:-${REPO_ROOT}/../results/chm13_chr22_index/chm13v2_chr22}"
GUIDE_TSV="${GUIDE_TSV:-${BENCH_DIR}/chr22_guides.tsv}"
PAM="${PAM:-GG}"
MAX_MISMATCHES="${MAX_MISMATCHES:-0}"
MAX_BULGES="${MAX_BULGES:-1}"
MAX_BULGE_SIZE="${MAX_BULGE_SIZE:-2}"
MIN_MATCH_FRACTION="${MIN_MATCH_FRACTION:-0.75}"
THREADS="${THREADS:-1}"
SAMPLE_INTERVAL="${SAMPLE_INTERVAL:-0.05}"
PYTHON="${PYTHON:-python3}"

CANDIDATE_E=$((MAX_MISMATCHES + MAX_BULGES * MAX_BULGE_SIZE))

if [[ -e "${RUN_DIR}" ]]; then
  echo "Refusing to overwrite existing run directory: ${RUN_DIR}" >&2
  exit 1
fi
mkdir -p "${RUN_DIR}/correctness" "${RUN_DIR}/timing" "${RUN_DIR}/tmp"

for path in "${CRISPRAPIDO_BIN}" "${COLUMBA_BIN}" "${REFERENCE}" "${GUIDE_TSV}"; do
  if [[ ! -e "${path}" ]]; then
    echo "Missing required path: ${path}" >&2
    exit 1
  fi
done
if [[ ! -x "${CRISPRAPIDO_BIN}" ]]; then
  echo "CRISPRapido binary is not executable: ${CRISPRAPIDO_BIN}" >&2
  exit 1
fi
if [[ ! -x "${COLUMBA_BIN}" ]]; then
  echo "Columba binary is not executable: ${COLUMBA_BIN}" >&2
  exit 1
fi
has_index=false
for suffix in meta bwt brt rev.brt sa.4 sa.bv.4 cct fsid sna pos headerSN.bin txt.bin; do
  if [[ -e "${COLUMBA_INDEX}.${suffix}" ]]; then
    has_index=true
    break
  fi
done
if [[ "${has_index}" != true ]]; then
  echo "Missing Columba index for prefix: ${COLUMBA_INDEX}" >&2
  exit 1
fi

sha_file() { sha256sum "$1" | awk '{print $1}'; }
count_lines() { if [[ -s "$1" ]]; then wc -l < "$1"; else echo 0; fi; }
metric_value() { awk -F'\t' -v key="$2" '$1 == key {print $2}' "$1"; }
index_size_bytes() { find "$(dirname "${COLUMBA_INDEX}")" -maxdepth 1 -type f -name "$(basename "${COLUMBA_INDEX}").*" -printf '%s\n' | awk '{s+=$1} END{print s+0}'; }
display_path() { realpath -m --relative-to="${REPO_ROOT}" "$1" 2>/dev/null || printf '%s' "$1"; }

write_environment() {
  {
    printf 'date_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'hostname\t%s\n' "$(hostname)"
    printf 'git_branch\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse --abbrev-ref HEAD)"
    printf 'git_commit\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD)"
    printf 'git_status_short\t%s\n' "$(git -C "${REPO_ROOT}" status --short | tr '\n' ';')"
    printf 'cpu_model\t%s\n' "$(awk -F': ' '/model name/ {print $2; exit}' /proc/cpuinfo 2>/dev/null || echo unknown)"
    printf 'logical_cpus\t%s\n' "$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo unknown)"
    printf 'total_ram_kib\t%s\n' "$(awk '/MemTotal/ {print $2}' /proc/meminfo 2>/dev/null || echo unknown)"
    printf 'rustc_version\t%s\n' "$(rustc +stable --version 2>/dev/null || rustc --version 2>/dev/null || echo unknown)"
    printf 'cargo_version\t%s\n' "$(cargo +stable --version 2>/dev/null || cargo --version 2>/dev/null || echo unknown)"
    printf 'crisprapido_bin\t%s\n' "$(display_path "${CRISPRAPIDO_BIN}")"
    printf 'crisprapido_sha256\t%s\n' "$(sha_file "${CRISPRAPIDO_BIN}")"
    printf 'crisprapido_mtime\t%s\n' "$(stat -c '%y' "${CRISPRAPIDO_BIN}")"
    printf 'columba_bin\t%s\n' "$(display_path "${COLUMBA_BIN}")"
    printf 'columba_sha256\t%s\n' "$(sha_file "${COLUMBA_BIN}")"
    printf 'columba_version\t%s\n' "$(${COLUMBA_BIN} --help 2>&1 | head -n 1)"
    printf 'reference\t%s\n' "$(display_path "${REFERENCE}")"
    printf 'reference_size_bytes\t%s\n' "$(stat -c '%s' "${REFERENCE}")"
    printf 'reference_sha256\t%s\n' "$(sha_file "${REFERENCE}")"
    printf 'columba_index_prefix\t%s\n' "$(display_path "${COLUMBA_INDEX}")"
    printf 'columba_index_size_bytes\t%s\n' "$(index_size_bytes)"
    printf 'guide_tsv\t%s\n' "$(display_path "${GUIDE_TSV}")"
    printf 'guide_tsv_sha256\t%s\n' "$(sha_file "${GUIDE_TSV}")"
    printf 'pam\t%s\n' "${PAM}"
    printf 'max_mismatches\t%s\n' "${MAX_MISMATCHES}"
    printf 'max_bulges\t%s\n' "${MAX_BULGES}"
    printf 'max_bulge_size\t%s\n' "${MAX_BULGE_SIZE}"
    printf 'candidate_e\t%s\n' "${CANDIDATE_E}"
    printf 'min_match_fraction\t%s\n' "${MIN_MATCH_FRACTION}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'sample_interval_seconds\t%s\n' "${SAMPLE_INTERVAL}"
    printf 'memory_method\t%s\n' "collect_peak_rss.py samples aggregate VmRSS for process plus descendants via /proc at fixed interval; very short-lived child peaks between samples may be missed. wait4/resource reports child user/system CPU seconds."
  } > "${RUN_DIR}/environment.txt"
}

make_command_file() {
  local display_crisprapido display_reference display_columba display_index
  display_crisprapido="$(display_path "${CRISPRAPIDO_BIN}")"
  display_reference="$(display_path "${REFERENCE}")"
  display_columba="$(display_path "${COLUMBA_BIN}")"
  display_index="$(display_path "${COLUMBA_INDEX}")"
  {
    printf 'mode\tcommand\n'
    printf 'baseline\t%q -r %q -g <GUIDE> -p %q -m %q -b %q -z %q -f %q -t %q\n' "${display_crisprapido}" "${display_reference}" "${PAM}" "${MAX_MISMATCHES}" "${MAX_BULGES}" "${MAX_BULGE_SIZE}" "${MIN_MATCH_FRACTION}" "${THREADS}"
    printf 'columba\t%q -r %q -g <GUIDE> -p %q -m %q -b %q -z %q -f %q -t %q --columba-bin %q --columba-index %q\n' "${display_crisprapido}" "${display_reference}" "${PAM}" "${MAX_MISMATCHES}" "${MAX_BULGES}" "${MAX_BULGE_SIZE}" "${MIN_MATCH_FRACTION}" "${THREADS}" "${display_columba}" "${display_index}"
  } > "${RUN_DIR}/commands.tsv"
}

run_one() {
  local phase="$1" mode="$2" guide_id="$3" guide_seq="$4" iter="$5"
  local out_dir="${RUN_DIR}/${phase}/${mode}/${iter}/${guide_id}"
  mkdir -p "${out_dir}"
  local paf="${out_dir}/output.paf"
  local stderr="${out_dir}/stderr.txt"
  local metrics="${out_dir}/metrics.tsv"
  local cmd=("${CRISPRAPIDO_BIN}" -r "${REFERENCE}" -g "${guide_seq}" -p "${PAM}" -m "${MAX_MISMATCHES}" -b "${MAX_BULGES}" -z "${MAX_BULGE_SIZE}" -f "${MIN_MATCH_FRACTION}" -t "${THREADS}")
  if [[ "${mode}" == "columba" ]]; then
    cmd+=(--columba-bin "${COLUMBA_BIN}" --columba-index "${COLUMBA_INDEX}")
  fi
  set +e
  env -u RUSTFLAGS -u LIBRARY_PATH -u LD_LIBRARY_PATH "${PYTHON}" "${BENCH_DIR}/collect_peak_rss.py" --stdout "${paf}" --stderr "${stderr}" --metrics "${metrics}" --sample-interval "${SAMPLE_INTERVAL}" -- "${cmd[@]}"
  local status=$?
  set -e
  printf '%q ' "${cmd[@]}" > "${out_dir}/command.txt"
  echo "${status}"
}

run_correctness_oracle() {
  RUN_DIR="${RUN_DIR}" REFERENCE="${REFERENCE}" GUIDE_TSV="${GUIDE_TSV}" PAM="${PAM}" "${PYTHON}" <<'PY'
import csv, os, re
from pathlib import Path
from collections import Counter
run=Path(os.environ['RUN_DIR']); ref_path=Path(os.environ['REFERENCE']); guide_tsv=Path(os.environ['GUIDE_TSV']); pam_pattern=os.environ['PAM'].upper()
comp=str.maketrans('ACGTNacgtn','TGCANtgcan')
def rc(s): return s.translate(comp)[::-1].upper()
def read_fasta(path):
    name=None; parts=[]
    for line in open(path):
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: break
            name=line[1:].split()[0]
        else: parts.append(line)
    return ''.join(parts).upper()
ref=read_fasta(ref_path)
guides=[]
with open(guide_tsv) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        guides.append((row['guide_id'], row['guide_sequence'].upper()))
def parse_paf(path,gid,guide):
    rows=[]
    if not path.exists(): return rows
    for line in open(path):
        if not line.strip(): continue
        f=line.rstrip('\n').split('\t'); tags={}
        for t in f[12:]:
            p=t.split(':',2)
            if len(p)==3: tags[p[0]]=p[2]
        rows.append({'guide_id':gid,'guide':guide,'rname':f[5],'strand':f[4],'start':int(f[7]),'end':int(f[8]),'cg':tags.get('cg',''),'nm':tags.get('nm',''),'ng':tags.get('ng',''),'bs':tags.get('bs',''),'cf':tags.get('cf',''),'line':line.rstrip('\n')})
    return rows
def exists_valid_alignment(guide,target):
    L=len(guide); T=len(target)
    if T==L and guide==target: return True,'20=','exact'
    if T in (L-1,L-2):
        k=L-T
        for p in range(0,L-k+1):
            if guide[:p]+guide[p+k:]==target:
                return True,(f'{p}=' if p else '')+f'{k}I'+(f'{L-p-k}=' if L-p-k else ''),'guide_insertion'
    if T in (L+1,L+2):
        k=T-L
        for p in range(0,L+1):
            if target[:p]==guide[:p] and target[p+k:]==guide[p:]:
                return True,(f'{p}=' if p else '')+f'{k}D'+(f'{L-p}=' if L-p else ''),'reference_deletion'
    return False,'','invalid'
def eval_row(r):
    raw=ref[r['start']:r['end']]
    target=rc(raw) if r['strand']=='-' else raw
    pam=ref[r['end']:r['end']+2] if r['strand']=='+' else (rc(ref[r['start']-2:r['start']]) if r['start']>=2 else '')
    ok,ocg,cls=exists_valid_alignment(r['guide'],target)
    return {**r,'target':target,'pam':pam,'oracle_valid':ok,'oracle_cigar':ocg,'oracle_class':cls,'pam_valid':pam==pam_pattern,'intended_valid':ok and pam==pam_pattern}
def equiv(a,b):
    if a['guide_id']!=b['guide_id'] or a['rname']!=b['rname'] or a['strand']!=b['strand'] or a['pam']!=b['pam']: return False
    if not a['intended_valid'] or not b['intended_valid']: return False
    if max(abs(a['start']-b['start']), abs(a['end']-b['end'])) > 2: return False
    return not (a['end'] <= b['start'] or b['end'] <= a['start'])
per=[]; all_b=[]; all_c=[]
for gid,guide in guides:
    b=[eval_row(r) for r in parse_paf(run/'correctness'/'baseline'/'pilot'/gid/'output.paf',gid,guide)]
    c=[eval_row(r) for r in parse_paf(run/'correctness'/'columba'/'pilot'/gid/'output.paf',gid,guide)]
    all_b += b; all_c += c
    valid_b=[x for x in b if x['intended_valid']]; valid_c=[x for x in c if x['intended_valid']]
    used=set(); shared=0
    for br in valid_b:
        hit=None
        for i,cr in enumerate(valid_c):
            if i not in used and equiv(br,cr): hit=i; break
        if hit is not None:
            used.add(hit); shared+=1
    per.append({'guide_id':gid,'baseline_exit':open(run/'correctness'/'baseline'/'pilot'/gid/'metrics.tsv').read().split('exit_status\t')[1].split('\n')[0],
                'columba_exit':open(run/'correctness'/'columba'/'pilot'/gid/'metrics.tsv').read().split('exit_status\t')[1].split('\n')[0],
                'baseline_raw_records':len(b),'columba_raw_records':len(c),'baseline_valid_loci':len(valid_b),'columba_valid_loci':len(valid_c),
                'shared_baseline_loci':shared,'baseline_missing_from_columba':len(valid_b)-shared,'columba_only_valid_loci':len(valid_c)-len(used),
                'columba_invalid_records':sum(not x['oracle_valid'] for x in c),'columba_non_gg_pam_records':sum(x['oracle_valid'] and not x['pam_valid'] for x in c),
                'columba_classes':dict(Counter(x['oracle_class'] for x in c))})
# aggregate by summing per guide, using guide-separated equivalence.
agg={k:sum(int(row[k]) for row in per) for k in ['baseline_raw_records','columba_raw_records','baseline_valid_loci','columba_valid_loci','shared_baseline_loci','baseline_missing_from_columba','columba_only_valid_loci','columba_invalid_records','columba_non_gg_pam_records']}
exit_ok=all(row['baseline_exit']=='0' and row['columba_exit']=='0' for row in per)
timing_eligible=exit_ok and agg['baseline_missing_from_columba']==0 and agg['columba_invalid_records']==0
with open(run/'per_guide_correctness.tsv','w',newline='') as fh:
    fields=['guide_id','baseline_exit','columba_exit','baseline_raw_records','columba_raw_records','baseline_valid_loci','columba_valid_loci','shared_baseline_loci','baseline_missing_from_columba','columba_only_valid_loci','columba_invalid_records','columba_non_gg_pam_records','columba_classes']
    w=csv.DictWriter(fh, fields, delimiter='\t'); w.writeheader(); w.writerows(per)
with open(run/'correctness_summary.tsv','w',newline='') as fh:
    fields=['baseline_raw_records','columba_raw_records','baseline_valid_loci','columba_valid_loci','shared_baseline_loci','baseline_missing_from_columba','columba_only_valid_loci','columba_invalid_records','columba_non_gg_pam_records','timing_eligible']
    w=csv.DictWriter(fh, fields, delimiter='\t'); w.writeheader(); w.writerow({**agg,'timing_eligible':'yes' if timing_eligible else 'no'})
# keep evaluated columba-only records for inspection without bloating top-level files.
with open(run/'correctness'/'oracle_columba_rows.tsv','w',newline='') as fh:
    fields=['guide_id','strand','start','end','cg','oracle_cigar','oracle_class','pam','pam_valid','oracle_valid','intended_valid','target']
    w=csv.DictWriter(fh, fields, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(all_c)
if not timing_eligible: raise SystemExit(2)
PY
}

write_timing_summary() {
  RUN_DIR="${RUN_DIR}" "${PYTHON}" <<'PY'
import csv, statistics, os
from pathlib import Path
run=Path(os.environ['RUN_DIR'])
rows=list(csv.DictReader(open(run/'run_summary.tsv'), delimiter='\t'))
measured=[r for r in rows if r['iteration_role']=='measured']
out=[]
for mode in ['baseline','columba']:
    vals=[r for r in measured if r['mode']==mode]
    walls=[float(r['total_wall_seconds']) for r in vals]
    rss=[int(r['max_peak_rss_kib']) for r in vals]
    recs=[int(r['total_paf_records']) for r in vals]
    out.append({'mode':mode,'measured_runs':len(vals),'wall_seconds':','.join(f'{x:.6f}' for x in walls),'median_wall_seconds':f'{statistics.median(walls):.6f}','min_wall_seconds':f'{min(walls):.6f}','max_wall_seconds':f'{max(walls):.6f}','peak_rss_kib':','.join(str(x) for x in rss),'median_peak_rss_kib':str(int(statistics.median(rss))),'total_paf_records':','.join(str(x) for x in recs),'deterministic_stdout':str(len(set(r['combined_stdout_sha256'] for r in vals))==1).lower()})
base=[x for x in out if x['mode']=='baseline'][0]; col=[x for x in out if x['mode']=='columba'][0]
speed=float(base['median_wall_seconds'])/float(col['median_wall_seconds'])
mem=float(col['median_peak_rss_kib'])/float(base['median_peak_rss_kib'])
with open(run/'timing_summary.tsv','w',newline='') as fh:
    fields=['mode','measured_runs','wall_seconds','median_wall_seconds','min_wall_seconds','max_wall_seconds','peak_rss_kib','median_peak_rss_kib','total_paf_records','deterministic_stdout']
    w=csv.DictWriter(fh, fields, delimiter='\t'); w.writeheader(); w.writerows(out)
    fh.write(f'speedup_baseline_median_over_columba\t{speed:.6f}\n')
    fh.write(f'memory_ratio_columba_over_baseline\t{mem:.6f}\n')
PY
}

printf 'guide_id\tmode\tphase\titeration\texit_status\twall_seconds\tuser_seconds\tsystem_seconds\tpeak_rss_kib\tpaf_records\tcandidate_count\tstdout_sha256\tstderr_sha256\n' > "${RUN_DIR}/guide_results.tsv"
printf 'run_label\tmode\titeration_role\titeration\tguide_count\texit_failures\ttotal_paf_records\ttotal_candidate_count\ttotal_wall_seconds\tmax_peak_rss_kib\tcombined_stdout_sha256\n' > "${RUN_DIR}/run_summary.tsv"

write_environment
make_command_file
cp "${RUN_DIR}/environment.txt" "${BENCH_DIR}/environment.txt"
cp "${RUN_DIR}/commands.tsv" "${BENCH_DIR}/commands.tsv"

mapfile -t guides < <(awk 'BEGIN{FS="\t"} NR>1 {print $1"\t"$2}' "${GUIDE_TSV}")

for row in "${guides[@]}"; do
  guide_id="${row%%$'\t'*}"
  guide_seq="${row#*$'\t'}"
  for mode in baseline columba; do
    status="$(run_one correctness "${mode}" "${guide_id}" "${guide_seq}" pilot)"
    m="${RUN_DIR}/correctness/${mode}/pilot/${guide_id}/metrics.tsv"
    paf="${RUN_DIR}/correctness/${mode}/pilot/${guide_id}/output.paf"
    printf '%s\t%s\tcorrectness\tpilot\t%s\t%s\t%s\t%s\t%s\t%s\tnot_exposed\t%s\t%s\n' \
      "${guide_id}" "${mode}" "${status}" \
      "$(metric_value "${m}" wall_seconds)" "$(metric_value "${m}" user_seconds)" "$(metric_value "${m}" system_seconds)" "$(metric_value "${m}" peak_rss_kib)" \
      "$(count_lines "${paf}")" "$(metric_value "${m}" stdout_sha256)" "$(metric_value "${m}" stderr_sha256)" >> "${RUN_DIR}/guide_results.tsv"
  done
done

if run_correctness_oracle; then
  correctness_ok=true
else
  correctness_ok=false
fi

cp "${RUN_DIR}/guide_results.tsv" "${BENCH_DIR}/guide_results.tsv"
cp "${RUN_DIR}/correctness_summary.tsv" "${BENCH_DIR}/correctness_summary.tsv"
cp "${RUN_DIR}/per_guide_correctness.tsv" "${BENCH_DIR}/per_guide_correctness.tsv"

run_panel() {
  local mode="$1" role="$2" iter="$3"
  local tmp_sha="${RUN_DIR}/timing/${mode}/${iter}/stdout.sha256.list"
  mkdir -p "$(dirname "${tmp_sha}")"
  : > "${tmp_sha}"
  local failures=0 total_records=0 total_wall=0 max_rss=0 count=0
  for row in "${guides[@]}"; do
    guide_id="${row%%$'\t'*}"
    guide_seq="${row#*$'\t'}"
    status="$(run_one timing "${mode}" "${guide_id}" "${guide_seq}" "${iter}")"
    m="${RUN_DIR}/timing/${mode}/${iter}/${guide_id}/metrics.tsv"
    paf="${RUN_DIR}/timing/${mode}/${iter}/${guide_id}/output.paf"
    records="$(count_lines "${paf}")"
    wall="$(metric_value "${m}" wall_seconds)"
    rss="$(metric_value "${m}" peak_rss_kib)"
    sha="$(metric_value "${m}" stdout_sha256)"
    if [[ "${status}" != "0" ]]; then failures=$((failures + 1)); fi
    total_records=$((total_records + records))
    total_wall="$(awk -v a="${total_wall}" -v b="${wall}" 'BEGIN{printf "%.6f", a+b}')"
    if [[ "${rss}" -gt "${max_rss}" ]]; then max_rss="${rss}"; fi
    echo "${guide_id} ${sha}" >> "${tmp_sha}"
    count=$((count + 1))
    printf '%s\t%s\ttiming\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tnot_exposed\t%s\t%s\n' \
      "${guide_id}" "${mode}" "${iter}" "${status}" "${wall}" "$(metric_value "${m}" user_seconds)" "$(metric_value "${m}" system_seconds)" "${rss}" "${records}" "${sha}" "$(metric_value "${m}" stderr_sha256)" >> "${RUN_DIR}/guide_results.tsv"
  done
  combined_sha="$(sha256sum "${tmp_sha}" | awk '{print $1}')"
  printf '%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\tnot_exposed\t%s\t%s\t%s\n' "${mode}" "${iter}" "${mode}" "${role}" "${iter}" "${count}" "${failures}" "${total_records}" "${total_wall}" "${max_rss}" "${combined_sha}" >> "${RUN_DIR}/run_summary.tsv"
}

if [[ "${correctness_ok}" != true ]]; then
  echo "Correctness eligibility failed; timing phase skipped. See ${RUN_DIR}/correctness_summary.tsv" >&2
  cp "${RUN_DIR}/run_summary.tsv" "${BENCH_DIR}/run_summary.tsv"
  exit 1
fi

run_panel baseline warmup warmup
run_panel columba warmup warmup
run_panel baseline measured measured_1
run_panel columba measured measured_1
run_panel columba measured measured_2
run_panel baseline measured measured_2
run_panel baseline measured measured_3
run_panel columba measured measured_3

write_timing_summary

cp "${RUN_DIR}/guide_results.tsv" "${BENCH_DIR}/guide_results.tsv"
cp "${RUN_DIR}/run_summary.tsv" "${BENCH_DIR}/run_summary.tsv"
cp "${RUN_DIR}/timing_summary.tsv" "${BENCH_DIR}/timing_summary.tsv"
cp "${RUN_DIR}/environment.txt" "${BENCH_DIR}/environment.txt"
cp "${RUN_DIR}/commands.tsv" "${BENCH_DIR}/commands.tsv"

echo "run_dir=${RUN_DIR}"
cat "${RUN_DIR}/correctness_summary.tsv"
echo "---"
cat "${RUN_DIR}/run_summary.tsv"
echo "---"
cat "${RUN_DIR}/timing_summary.tsv"
