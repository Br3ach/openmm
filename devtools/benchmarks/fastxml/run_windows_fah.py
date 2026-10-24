#!/usr/bin/env python3
"""Validate copied FAH XML and compare native Windows serializer read measurements."""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import re
import statistics
import subprocess
import sys


def file_hash(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def normalized_hash(path, version_attribute):
    with path.open("rb") as stream:
        prefix = stream.read(4096)
        attributes = list(re.finditer(rb'openmmVersion="[^"]*"', prefix))
        if len(attributes) != 1:
            raise RuntimeError("Expected one root OpenMM version attribute")
        match = attributes[0]
        digest = hashlib.sha256(prefix[:match.start()] + version_attribute + prefix[match.end():])
        while chunk := stream.read(4*1024*1024):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--system", type=Path)
    parser.add_argument("--state", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--memory-limit-mib", type=int, default=6144)
    args = parser.parse_args()
    if sys.platform != "win32":
        parser.error("This harness requires native Windows")
    if args.repeats < 1 or not 1 <= args.memory_limit_mib <= 65536:
        parser.error("Invalid repeat count or memory limit")
    cases = [(kind, path.resolve()) for kind, path in [("system", args.system), ("state", args.state)] if path]
    if not cases:
        parser.error("Provide --system and/or --state using stable, uncompressed copies")
    binaries = {"candidate": args.candidate.resolve()}
    if args.baseline:
        binaries["baseline"] = args.baseline.resolve()
    args.output.mkdir(parents=True, exist_ok=True)

    def execute(variant, kind, source, target):
        result = subprocess.run([str(binaries[variant]), kind, str(source), str(target),
                                 "1", str(args.memory_limit_mib)],
                                capture_output=True, text=True, timeout=180)
        if result.returncode:
            raise RuntimeError(f"{variant} {kind}: {result.returncode}: {result.stderr}")
        return json.loads(result.stdout)

    manifests = {}
    verification = []
    groups = collections.defaultdict(list)
    with (args.output/"results.jsonl").open("w") as output:
        for kind, source in cases:
            original_stat = source.stat()
            manifests[kind] = {"bytes": original_stat.st_size, "sha256": file_hash(source)}
            # Each run reads the same file, with its pages warmed outside the measured interval.
            for repeat in range(args.repeats):
                order = list(binaries)
                if repeat % 2:
                    order.reverse()
                for variant in order:
                    with source.open("rb") as data:
                        while data.read(4*1024*1024):
                            pass
                    row = execute(variant, kind, source, "-")
                    row.update(variant=variant, repeat=repeat+1)
                    groups[(kind, variant)].append(row)
                    output.write(json.dumps(row)+"\n")
                    output.flush()
                    print(json.dumps(row), flush=True)
            # The complete output fingerprint checks every serialized value and force.
            target = args.output/(kind+"-roundtrip.xml")
            if target.resolve() == source:
                raise RuntimeError("Round-trip output would overwrite the input")
            row = execute("candidate", kind, source, target)
            with target.open("rb") as data:
                prefix = data.read(4096)
            attribute = re.search(rb'openmmVersion="[^"]*"', prefix)
            if not attribute:
                raise RuntimeError("Output has no root OpenMM version attribute")
            reference = normalized_hash(source, attribute.group())
            actual = file_hash(target)
            identical = reference == actual
            verification.append(dict(row, input_bytes=source.stat().st_size,
                                     output_bytes=target.stat().st_size,
                                     output_sha256=actual, normalized_input_sha256=reference,
                                     equal_after_root_version_normalization=identical))
            (args.output/"verification.json").write_text(json.dumps(verification, indent=2)+"\n")
            if not identical:
                raise RuntimeError("XML differs beyond its root OpenMM version attribute")
            if (original_stat.st_size, original_stat.st_mtime_ns) != (source.stat().st_size, source.stat().st_mtime_ns):
                raise RuntimeError("Input changed during validation")
            target.unlink()
    summary = []
    for (kind, variant), rows in groups.items():
        counts = {row["particles"] for row in rows}
        if len(counts) != 1:
            raise RuntimeError("Particle count changed across repeated runs")
        summary.append(dict(kind=kind, variant=variant, particles=rows[0]["particles"],
                            samples=len(rows), read_seconds_median=statistics.median(r["read_seconds"] for r in rows),
                            peak_working_set_kib_median=statistics.median(r["read_peak_working_set_kib"] for r in rows),
                            peak_scope="process_lifetime"))
    for kind, _ in cases:
        counts = {row["particles"] for (case, _), rows in groups.items() if case == kind for row in rows}
        if len(counts) != 1:
            raise RuntimeError("Particle count differs between variants")
    (args.output/"fixtures.json").write_text(json.dumps(manifests, indent=2)+"\n")
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")


if __name__ == "__main__":
    main()
