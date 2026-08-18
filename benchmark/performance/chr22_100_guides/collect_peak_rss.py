#!/usr/bin/env python3
import argparse
import hashlib
import os
import resource
import subprocess
import sys
import time
from pathlib import Path


def read_status_kib(pid: int, key: str):
    try:
        with open(f"/proc/{pid}/status", "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(key + ":"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return int(parts[1])
    except FileNotFoundError:
        return None
    except ProcessLookupError:
        return None
    return None


def descendants(pid: int):
    children = []
    task_children = Path(f"/proc/{pid}/task")
    try:
        tasks = list(task_children.iterdir())
    except FileNotFoundError:
        return children
    for task in tasks:
        try:
            data = (task / "children").read_text().strip()
        except FileNotFoundError:
            continue
        if not data:
            continue
        for child in data.split():
            try:
                child_pid = int(child)
            except ValueError:
                continue
            children.append(child_pid)
            children.extend(descendants(child_pid))
    return children


def aggregate_rss_kib(pid: int):
    total = 0
    seen = set()
    for proc_pid in [pid] + descendants(pid):
        if proc_pid in seen:
            continue
        seen.add(proc_pid)
        rss = read_status_kib(proc_pid, "VmRSS")
        if rss is not None:
            total += rss
    return total


def sha256_file(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Run a command and sample peak aggregate RSS from /proc.")
    parser.add_argument("--stdout", required=True)
    parser.add_argument("--stderr", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--sample-interval", type=float, default=0.05)
    parser.add_argument("--", dest="separator", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise SystemExit("missing command")

    stdout_path = Path(args.stdout)
    stderr_path = Path(args.stderr)
    metrics_path = Path(args.metrics)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)

    start_wall = time.monotonic()
    start_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_rss_kib = 0

    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        proc = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        while True:
            rc = proc.poll()
            rss = aggregate_rss_kib(proc.pid)
            if rss > peak_rss_kib:
                peak_rss_kib = rss
            if rc is not None:
                break
            time.sleep(args.sample_interval)
        exit_status = proc.returncode

    end_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    elapsed = time.monotonic() - start_wall
    user_seconds = end_usage.ru_utime - start_usage.ru_utime
    sys_seconds = end_usage.ru_stime - start_usage.ru_stime

    metrics = {
        "exit_status": str(exit_status),
        "wall_seconds": f"{elapsed:.6f}",
        "user_seconds": f"{user_seconds:.6f}",
        "system_seconds": f"{sys_seconds:.6f}",
        "peak_rss_kib": str(peak_rss_kib),
        "stdout_bytes": str(stdout_path.stat().st_size if stdout_path.exists() else 0),
        "stderr_bytes": str(stderr_path.stat().st_size if stderr_path.exists() else 0),
        "stdout_sha256": sha256_file(stdout_path) if stdout_path.exists() else "missing",
        "stderr_sha256": sha256_file(stderr_path) if stderr_path.exists() else "missing",
    }
    with metrics_path.open("w", encoding="utf-8") as handle:
        for key, value in metrics.items():
            handle.write(f"{key}\t{value}\n")

    return exit_status


if __name__ == "__main__":
    sys.exit(main())
