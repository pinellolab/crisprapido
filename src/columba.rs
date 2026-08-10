use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Debug)]
pub(crate) struct ColumbaRunConfig<'a> {
    pub(crate) columba_bin: &'a Path,
    pub(crate) index_prefix: &'a Path,
    pub(crate) guide: &'a str,
    pub(crate) candidate_edit_distance: u32,
    pub(crate) threads: Option<usize>,
    pub(crate) keep_sam: bool,
}

#[derive(Debug)]
pub(crate) struct ColumbaRunOutput {
    pub(crate) sam_path: PathBuf,
    temp_dir: PathBuf,
    keep_files: bool,
    cleaned: bool,
}

impl ColumbaRunOutput {
    #[cfg(test)]
    pub(crate) fn temp_dir(&self) -> &Path {
        &self.temp_dir
    }

    fn cleanup(&mut self) {
        if self.cleaned || self.keep_files {
            return;
        }
        if let Err(e) = fs::remove_dir_all(&self.temp_dir) {
            eprintln!(
                "Warning: failed to remove temporary Columba directory '{}': {}",
                self.temp_dir.display(),
                e
            );
        }
        self.cleaned = true;
    }
}

impl Drop for ColumbaRunOutput {
    fn drop(&mut self) {
        self.cleanup();
    }
}

fn temp_columba_dir() -> PathBuf {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        "crisprapido-columba-{}-{}",
        std::process::id(),
        now
    ))
}

fn prefixed_index_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value: OsString = prefix.as_os_str().to_os_string();
    value.push(".");
    value.push(suffix);
    PathBuf::from(value)
}

fn columba_index_exists(prefix: &Path) -> bool {
    if prefix.is_file() {
        return true;
    }

    [
        "meta",
        "bwt",
        "brt",
        "rev.brt",
        "sa.4",
        "sa.bv.4",
        "cct",
        "fsid",
        "sna",
        "pos",
        "headerSN.bin",
        "txt.bin",
    ]
    .iter()
    .any(|suffix| prefixed_index_path(prefix, suffix).exists())
}

#[cfg(unix)]
fn is_executable(path: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;

    path.metadata()
        .map(|metadata| metadata.is_file() && metadata.permissions().mode() & 0o111 != 0)
        .unwrap_or(false)
}

#[cfg(not(unix))]
fn is_executable(path: &Path) -> bool {
    path.is_file()
}

fn write_columba_log(path: &Path, stdout: &[u8], stderr: &[u8]) -> Result<(), String> {
    let mut file = File::create(path)
        .map_err(|e| format!("Failed to create Columba log '{}': {}", path.display(), e))?;
    file.write_all(b"stdout:\n")
        .and_then(|_| file.write_all(stdout))
        .and_then(|_| file.write_all(b"\nstderr:\n"))
        .and_then(|_| file.write_all(stderr))
        .map_err(|e| format!("Failed to write Columba log '{}': {}", path.display(), e))
}

pub(crate) fn candidate_edit_distance_bound(
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
) -> u32 {
    // WFA2 remains the authoritative filter: this total edit bound is a
    // candidate superset for at most m substitutions plus b gap groups of
    // at most z bases each.
    max_mismatches.saturating_add(max_bulges.saturating_mul(max_bulge_size))
}

pub(crate) fn run_columba(config: &ColumbaRunConfig<'_>) -> Result<ColumbaRunOutput, String> {
    if !config.columba_bin.exists() {
        return Err(format!(
            "Columba executable not found: {}",
            config.columba_bin.display()
        ));
    }
    if !is_executable(config.columba_bin) {
        return Err(format!(
            "Columba path is not executable: {}",
            config.columba_bin.display()
        ));
    }
    if !columba_index_exists(config.index_prefix) {
        return Err(format!(
            "Columba index not found for prefix: {}",
            config.index_prefix.display()
        ));
    }

    let temp_dir = temp_columba_dir();
    fs::create_dir(&temp_dir).map_err(|e| {
        format!(
            "Failed to create temporary Columba directory '{}': {}",
            temp_dir.display(),
            e
        )
    })?;
    let query_path = temp_dir.join("guide.fa");
    let sam_path = temp_dir.join("output.sam");
    let log_path = temp_dir.join("columba.log");

    let run_result = (|| {
        fs::write(&query_path, format!(">guide\n{}\n", config.guide)).map_err(|e| {
            format!(
                "Failed to write temporary Columba query FASTA '{}': {}",
                query_path.display(),
                e
            )
        })?;

        let mut command = Command::new(config.columba_bin);
        command
            .arg("-r")
            .arg(config.index_prefix)
            .arg("-f")
            .arg(&query_path)
            .arg("-a")
            .arg("all")
            .arg("-m")
            .arg("edit")
            .arg("-e")
            .arg(config.candidate_edit_distance.to_string());
        if let Some(threads) = config.threads {
            command.arg("-t").arg(threads.to_string());
        }
        command.arg("-o").arg(&sam_path);

        let output = command.output().map_err(|e| {
            format!(
                "Failed to execute Columba '{}': {}",
                config.columba_bin.display(),
                e
            )
        })?;
        write_columba_log(&log_path, &output.stdout, &output.stderr)?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            if !stderr.trim().is_empty() {
                eprintln!("{}", stderr.trim_end());
            }
            return Err(format!(
                "Columba exited with status {}. See log: {}",
                output.status,
                log_path.display()
            ));
        }

        if !sam_path.exists() {
            return Err(format!(
                "Columba completed successfully but did not create SAM output: {}",
                sam_path.display()
            ));
        }

        Ok(ColumbaRunOutput {
            sam_path,
            temp_dir: temp_dir.clone(),
            keep_files: config.keep_sam,
            cleaned: false,
        })
    })();

    if run_result.is_err() && !config.keep_sam {
        let _ = fs::remove_dir_all(&temp_dir);
    }

    run_result
}

#[derive(Debug)]
pub(crate) struct ColumbaCandidate {
    pub(crate) query_name: String,
    pub(crate) reference_name: String,
    pub(crate) reference_start: usize,
    pub(crate) reverse: bool,
    pub(crate) cigar: String,
    pub(crate) edit_distance: usize,
    pub(crate) alignment_score: i32,
}

pub(crate) fn deduplicate_columba_candidates(
    candidates: Vec<ColumbaCandidate>,
) -> Vec<ColumbaCandidate> {
    let mut seen = HashSet::new();
    let mut deduplicated = Vec::new();

    for candidate in candidates {
        let key = (
            candidate.query_name.clone(),
            candidate.reference_name.clone(),
            candidate.reference_start,
            candidate.reverse,
            candidate.cigar.clone(),
        );
        if seen.insert(key) {
            deduplicated.push(candidate);
        }
    }

    deduplicated
}

pub(crate) fn run_columba_candidate_generation(
    config: &ColumbaRunConfig<'_>,
) -> Result<Vec<ColumbaCandidate>, String> {
    let columba_output = run_columba(config)?;
    let mut candidates = parse_columba_sam_file(&columba_output.sam_path)
        .map_err(|e| format!("Failed to parse generated Columba SAM file: {}", e))?;

    if config.candidate_edit_distance > 0 {
        let exact_config = ColumbaRunConfig {
            candidate_edit_distance: 0,
            ..*config
        };
        let exact_output = run_columba(&exact_config)?;
        candidates.extend(
            parse_columba_sam_file(&exact_output.sam_path)
                .map_err(|e| format!("Failed to parse generated exact Columba SAM file: {}", e))?,
        );
        candidates = deduplicate_columba_candidates(candidates);
    }

    Ok(candidates)
}

pub(crate) fn cigar_reference_span(cigar: &str) -> Result<usize, String> {
    let mut span = 0;
    let mut length = 0usize;
    let mut has_op = false;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add(c.to_digit(10).unwrap() as usize))
                .ok_or_else(|| format!("CIGAR operation length overflow in '{}'", cigar))?;
            continue;
        }

        if length == 0 {
            return Err(format!(
                "Invalid CIGAR '{}': missing operation length",
                cigar
            ));
        }

        match c {
            'M' | 'D' | 'N' | '=' | 'X' => span += length,
            'I' | 'S' | 'H' | 'P' => {}
            _ => {
                return Err(format!(
                    "Invalid CIGAR '{}': unsupported operation '{}'",
                    cigar, c
                ))
            }
        }

        has_op = true;
        length = 0;
    }

    if length != 0 {
        return Err(format!(
            "Invalid CIGAR '{}': trailing operation length",
            cigar
        ));
    }

    if !has_op {
        return Err(format!("Invalid CIGAR '{}': no operations", cigar));
    }

    Ok(span)
}

pub(crate) fn parse_columba_sam_record(line: &str) -> Result<Option<ColumbaCandidate>, String> {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with('@') {
        return Ok(None);
    }

    let fields: Vec<&str> = trimmed.split('\t').collect();
    if fields.len() < 11 {
        return Err(format!(
            "Invalid SAM record: expected at least 11 columns, found {}",
            fields.len()
        ));
    }

    let flag = fields[1]
        .parse::<u16>()
        .map_err(|e| format!("Invalid SAM FLAG '{}': {}", fields[1], e))?;
    if flag & 4 != 0 {
        return Ok(None);
    }

    let pos = fields[3]
        .parse::<usize>()
        .map_err(|e| format!("Invalid SAM POS '{}': {}", fields[3], e))?;
    if pos == 0 {
        return Err("Invalid SAM POS '0': mapped records must be 1-based".to_string());
    }

    fields[4]
        .parse::<u16>()
        .map_err(|e| format!("Invalid SAM MAPQ '{}': {}", fields[4], e))?;
    let _reference_span = cigar_reference_span(fields[5])?;

    let mut edit_distance = None;
    let mut alignment_score = None;
    for tag in &fields[11..] {
        if let Some(value) = tag.strip_prefix("NM:i:") {
            edit_distance = Some(
                value
                    .parse::<usize>()
                    .map_err(|e| format!("Invalid NM tag '{}': {}", tag, e))?,
            );
        } else if let Some(value) = tag.strip_prefix("AS:i:") {
            alignment_score = Some(
                value
                    .parse::<i32>()
                    .map_err(|e| format!("Invalid AS tag '{}': {}", tag, e))?,
            );
        }
    }

    let candidate = ColumbaCandidate {
        query_name: fields[0].to_string(),
        reference_name: fields[2].to_string(),
        reference_start: pos - 1,
        reverse: flag & 16 != 0,
        cigar: fields[5].to_string(),
        edit_distance: edit_distance.ok_or_else(|| {
            format!(
                "Missing NM:i tag for mapped SAM record '{}' at {}:{}",
                fields[0], fields[2], fields[3]
            )
        })?,
        alignment_score: alignment_score.ok_or_else(|| {
            format!(
                "Missing AS:i tag for mapped SAM record '{}' at {}:{}",
                fields[0], fields[2], fields[3]
            )
        })?,
    };

    Ok(Some(candidate))
}

pub(crate) fn parse_columba_sam_file(path: &PathBuf) -> Result<Vec<ColumbaCandidate>, String> {
    let file = File::open(path)
        .map_err(|e| format!("Failed to open Columba SAM '{}': {}", path.display(), e))?;
    let reader = BufReader::new(file);
    let mut candidates = Vec::new();

    for (line_idx, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| {
            format!(
                "Failed to read Columba SAM '{}' line {}: {}",
                path.display(),
                line_idx + 1,
                e
            )
        })?;
        if let Some(candidate) = parse_columba_sam_record(&line)
            .map_err(|e| format!("{} on line {}", e, line_idx + 1))?
        {
            candidates.push(candidate);
        }
    }

    Ok(candidates)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_candidate(line: &str) -> ColumbaCandidate {
        parse_columba_sam_record(line).unwrap().unwrap()
    }

    #[test]
    fn test_parse_columba_primary_forward_record() {
        let candidate = parse_candidate(
            "guide_20bp\t0\tperfect_forward\t11\t60\t20M\t*\t0\t0\tGAGTCCGAGCAGAAGAAGAA\t*\tAS:i:0\tNM:i:0\tPG:Z:Columba",
        );

        assert_eq!(candidate.query_name, "guide_20bp");
        assert_eq!(candidate.reference_name, "perfect_forward");
        assert_eq!(candidate.reference_start, 10);
        assert!(!candidate.reverse);
        assert_eq!(candidate.cigar, "20M");
        assert_eq!(candidate.edit_distance, 0);
        assert_eq!(candidate.alignment_score, 0);
    }

    #[test]
    fn test_parse_columba_secondary_forward_record() {
        let candidate = parse_candidate(
            "guide_20bp\t256\tone_substitution\t11\t0\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1\tPG:Z:Columba",
        );

        assert_eq!(candidate.reference_name, "one_substitution");
        assert_eq!(candidate.reference_start, 10);
        assert!(!candidate.reverse);
        assert_eq!(candidate.edit_distance, 1);
        assert_eq!(candidate.alignment_score, 1);
    }

    #[test]
    fn test_parse_columba_secondary_reverse_record() {
        let candidate = parse_candidate(
            "guide_20bp\t272\tperfect_reverse_complement\t11\t0\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0\tPG:Z:Columba",
        );

        assert_eq!(candidate.reference_name, "perfect_reverse_complement");
        assert_eq!(candidate.reference_start, 10);
        assert!(candidate.reverse);
    }

    #[test]
    fn test_columba_pos_11_converts_to_start0_10() {
        let candidate = parse_candidate(
            "guide_20bp\t0\tchr1\t11\t60\t20M\t*\t0\t0\tGAGTCCGAGCAGAAGAAGAA\t*\tAS:i:0\tNM:i:0",
        );

        assert_eq!(candidate.reference_start, 10);
    }

    #[test]
    fn test_columba_cigar_reference_spans() {
        assert_eq!(cigar_reference_span("20M").unwrap(), 20);
        assert_eq!(cigar_reference_span("10M1D10M").unwrap(), 21);
        assert_eq!(cigar_reference_span("10M1I9M").unwrap(), 19);
    }

    #[test]
    fn test_parse_columba_ignores_unmapped_record() {
        let candidate =
            parse_columba_sam_record("guide_20bp\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0")
                .unwrap();

        assert!(candidate.is_none());
    }

    #[test]
    fn test_parse_columba_accepts_seq_star() {
        let candidate =
            parse_candidate("guide_20bp\t256\tchr1\t11\t0\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        assert_eq!(candidate.query_name, "guide_20bp");
        assert_eq!(candidate.cigar, "20M");
    }

    #[test]
    fn test_parse_columba_missing_nm_is_clear_error() {
        let error = parse_columba_sam_record(
            "guide_20bp\t0\tchr1\t11\t60\t20M\t*\t0\t0\tGAGTCCGAGCAGAAGAAGAA\t*\tAS:i:0",
        )
        .unwrap_err();

        assert!(error.contains("Missing NM:i tag"));
        assert!(error.contains("guide_20bp"));
        assert!(error.contains("chr1:11"));
    }
    fn unique_test_dir(name: &str) -> PathBuf {
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|duration| duration.as_nanos())
            .unwrap_or(0);
        let path = std::env::temp_dir().join(format!(
            "crisprapido-columba-test-{}-{}-{}",
            name,
            std::process::id(),
            now
        ));
        fs::create_dir(&path).unwrap();
        path
    }

    fn write_dummy_index(prefix: &Path) {
        fs::write(prefixed_index_path(prefix, "meta"), b"index").unwrap();
    }

    #[cfg(unix)]
    fn write_mock_executable(path: &Path, body: &str) {
        use std::os::unix::fs::PermissionsExt;

        fs::write(path, body).unwrap();
        let mut permissions = fs::metadata(path).unwrap().permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(path, permissions).unwrap();
    }

    #[test]
    fn test_candidate_edit_distance_bound_zero() {
        assert_eq!(candidate_edit_distance_bound(0, 0, 0), 0);
    }

    #[test]
    fn test_candidate_edit_distance_bound_single_two_base_bulge() {
        assert_eq!(candidate_edit_distance_bound(0, 1, 2), 2);
    }

    #[test]
    fn test_candidate_edit_distance_bound_mismatches_plus_bulge() {
        assert_eq!(candidate_edit_distance_bound(2, 1, 2), 4);
    }

    #[test]
    fn test_candidate_edit_distance_bound_multiple_gap_groups() {
        assert_eq!(candidate_edit_distance_bound(1, 2, 3), 7);
    }

    #[test]
    fn test_candidate_edit_distance_bound_saturates_on_overflow() {
        assert_eq!(candidate_edit_distance_bound(u32::MAX - 1, 2, 3), u32::MAX);
    }

    #[cfg(unix)]
    #[test]
    fn test_run_columba_missing_executable() {
        let dir = unique_test_dir("missing-executable");
        let index_prefix = dir.join("idx");
        write_dummy_index(&index_prefix);

        let error = run_columba(&ColumbaRunConfig {
            columba_bin: &dir.join("missing-columba"),
            index_prefix: &index_prefix,
            guide: "GAGTCCGAGCAGAAGAAGAA",
            candidate_edit_distance: 1,
            threads: None,
            keep_sam: false,
        })
        .unwrap_err();

        assert!(error.contains("Columba executable not found"));
        fs::remove_dir_all(dir).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn test_run_columba_missing_index() {
        let dir = unique_test_dir("missing-index");
        let bin = dir.join("columba");
        write_mock_executable(&bin, "#!/bin/sh\nexit 0\n");

        let error = run_columba(&ColumbaRunConfig {
            columba_bin: &bin,
            index_prefix: &dir.join("missing-index"),
            guide: "GAGTCCGAGCAGAAGAAGAA",
            candidate_edit_distance: 1,
            threads: None,
            keep_sam: false,
        })
        .unwrap_err();

        assert!(error.contains("Columba index not found"));
        fs::remove_dir_all(dir).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn test_run_columba_failure_reports_status_and_cleans_up() {
        let dir = unique_test_dir("failure");
        let bin = dir.join("columba");
        let index_prefix = dir.join("idx");
        write_dummy_index(&index_prefix);
        write_mock_executable(&bin, "#!/bin/sh\necho simulated failure >&2\nexit 7\n");

        let error = run_columba(&ColumbaRunConfig {
            columba_bin: &bin,
            index_prefix: &index_prefix,
            guide: "GAGTCCGAGCAGAAGAAGAA",
            candidate_edit_distance: 1,
            threads: Some(2),
            keep_sam: false,
        })
        .unwrap_err();

        assert!(error.contains("Columba exited with status"));
        fs::remove_dir_all(dir).unwrap();
    }

    #[cfg(unix)]
    fn successful_mock_script() -> &'static str {
        "#!/bin/sh\nout=\"\"\nwhile [ \"$#\" -gt 0 ]; do\n  if [ \"$1\" = \"-o\" ]; then\n    shift\n    out=\"$1\"\n  fi\n  shift\ndone\nprintf '@HD\\tVN:1.6\\n' > \"$out\"\nprintf 'guide\\t0\\tchr1\\t11\\t60\\t20M\\t*\\t0\\t0\\t*\\t*\\tAS:i:0\\tNM:i:0\\n' >> \"$out\"\necho columba stdout\necho columba stderr >&2\nexit 0\n"
    }

    #[cfg(unix)]
    #[test]
    fn test_run_columba_uses_candidate_edit_distance() {
        let dir = unique_test_dir("candidate-edit-distance");
        let bin = dir.join("columba");
        let index_prefix = dir.join("idx");
        write_dummy_index(&index_prefix);
        write_mock_executable(
            &bin,
            r#"#!/bin/sh
out=""
e=""
while [ "$#" -gt 0 ]; do
  if [ "$1" = "-e" ]; then
    shift
    e="$1"
  elif [ "$1" = "-o" ]; then
    shift
    out="$1"
  fi
  shift
done
printf '@HD	VN:1.6
' > "$out"
if [ "$e" = "0" ]; then
  printf 'guide	0	chr1	101	60	20M	*	0	0	*	*	AS:i:0	NM:i:0
' >> "$out"
else
  printf 'guide	0	chr1	99	60	20M	*	0	0	*	*	AS:i:1	NM:i:1
' >> "$out"
fi
exit 0
"#,
        );

        let candidates = run_columba_candidate_generation(&ColumbaRunConfig {
            columba_bin: &bin,
            index_prefix: &index_prefix,
            guide: "GAGTCCGAGCAGAAGAAGAA",
            candidate_edit_distance: 2,
            threads: None,
            keep_sam: false,
        })
        .unwrap();

        assert_eq!(candidates.len(), 2);
        assert!(candidates
            .iter()
            .any(|candidate| candidate.reference_start == 100 && candidate.edit_distance == 0));
        assert!(candidates
            .iter()
            .any(|candidate| candidate.reference_start == 98 && candidate.edit_distance == 1));
        fs::remove_dir_all(dir).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn test_run_columba_success_with_mocked_executable() {
        let dir = unique_test_dir("success");
        let bin = dir.join("columba");
        let index_prefix = dir.join("idx");
        write_dummy_index(&index_prefix);
        write_mock_executable(&bin, successful_mock_script());

        let output = run_columba(&ColumbaRunConfig {
            columba_bin: &bin,
            index_prefix: &index_prefix,
            guide: "GAGTCCGAGCAGAAGAAGAA",
            candidate_edit_distance: 1,
            threads: Some(2),
            keep_sam: false,
        })
        .unwrap();

        assert!(output.sam_path.exists());
        let candidates = parse_columba_sam_file(&output.sam_path).unwrap();
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].reference_name, "chr1");
        drop(output);
        fs::remove_dir_all(dir).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn test_run_columba_cleanup_removes_temp_files() {
        let dir = unique_test_dir("cleanup");
        let bin = dir.join("columba");
        let index_prefix = dir.join("idx");
        write_dummy_index(&index_prefix);
        write_mock_executable(&bin, successful_mock_script());

        let output = run_columba(&ColumbaRunConfig {
            columba_bin: &bin,
            index_prefix: &index_prefix,
            guide: "GAGTCCGAGCAGAAGAAGAA",
            candidate_edit_distance: 1,
            threads: None,
            keep_sam: false,
        })
        .unwrap();
        let temp_dir = output.temp_dir().to_path_buf();
        assert!(temp_dir.exists());
        drop(output);

        assert!(!temp_dir.exists());
        fs::remove_dir_all(dir).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn test_run_columba_keep_sam_preserves_temp_files() {
        let dir = unique_test_dir("keep");
        let bin = dir.join("columba");
        let index_prefix = dir.join("idx");
        write_dummy_index(&index_prefix);
        write_mock_executable(&bin, successful_mock_script());

        let output = run_columba(&ColumbaRunConfig {
            columba_bin: &bin,
            index_prefix: &index_prefix,
            guide: "GAGTCCGAGCAGAAGAAGAA",
            candidate_edit_distance: 1,
            threads: None,
            keep_sam: true,
        })
        .unwrap();
        let temp_dir = output.temp_dir().to_path_buf();
        let sam_path = output.sam_path.clone();
        drop(output);

        assert!(temp_dir.exists());
        assert!(sam_path.exists());
        fs::remove_dir_all(temp_dir).unwrap();
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn test_deduplicate_columba_candidates_retains_exact_e0_candidate() {
        let exact = parse_candidate("guide	0	chr1	101	60	20M	*	0	0	*	*	AS:i:0	NM:i:0");
        let shifted = parse_candidate("guide	0	chr1	99	60	20M	*	0	0	*	*	AS:i:1	NM:i:1");
        let duplicate_exact = parse_candidate("guide	256	chr1	101	60	20M	*	0	0	*	*	AS:i:0	NM:i:0");

        let candidates = deduplicate_columba_candidates(vec![shifted, exact, duplicate_exact]);

        assert_eq!(candidates.len(), 2);
        assert!(candidates
            .iter()
            .any(|candidate| candidate.reference_start == 100 && candidate.edit_distance == 0));
        assert!(candidates
            .iter()
            .any(|candidate| candidate.reference_start == 98 && candidate.edit_distance == 1));
    }
}
