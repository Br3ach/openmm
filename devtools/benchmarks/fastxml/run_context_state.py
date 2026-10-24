#!/usr/bin/env python3
"""Compare native Windows Context::getState with one caller and two DLL builds."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if sys.platform != 'win32' or args.repeats < 1:
        parser.error('Requires native Windows and at least one repeat')
    binaries = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
    fingerprints = {hashlib.sha256(p.read_bytes()).hexdigest() for p in binaries.values()}
    if len(fingerprints) != 1:
        parser.error('Use identical benchmark executables beside the two different OpenMM DLLs')
    args.output.mkdir(parents=True, exist_ok=True)
    cases = [(1000000, 3), (1000000, 31), (7950256, 3), (7950256, 31)]
    rows = []
    for count, mask in cases:
        for repeat in range(args.repeats):
            order = ['baseline', 'candidate'] if repeat % 2 == 0 else ['candidate', 'baseline']
            for variant in order:
                exe = binaries[variant]
                result = subprocess.run([str(exe), str(count), str(mask)], cwd=exe.parent,
                                        capture_output=True, text=True, timeout=120)
                if result.returncode:
                    raise RuntimeError(f'{variant}: {result.returncode}: {result.stderr}')
                row = json.loads(result.stdout)
                row.update(variant=variant, repeat=repeat+1)
                rows.append(row)
                (args.output/'raw.json').write_text(json.dumps(rows, indent=2)+'\n')
                print(json.dumps(row), flush=True)
    summary = []
    for count, mask in cases:
        selected = [r for r in rows if r['particles'] == count and r['mask'] == mask]
        if len({r['fingerprint'] for r in selected}) != 1:
            raise RuntimeError('Returned particle data or energy differs between DLLs')
        item = dict(particles=count, mask=mask, identical_fingerprints=True)
        for variant in ['baseline', 'candidate']:
            samples = [r for r in selected if r['variant'] == variant]
            item[variant] = {key: statistics.median(r[key] for r in samples)
                             for key in ['seconds', 'peak_working_set_kib', 'baseline_working_set_kib',
                                         'retained_working_set_kib', 'private_kib']}
        summary.append(item)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
