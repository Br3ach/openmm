#!/usr/bin/env python3
"""Compare Context construction using one Windows caller and two DLLs."""
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
    if len({hashlib.sha256(p.read_bytes()).hexdigest() for p in binaries.values()}) != 1:
        parser.error('Use identical caller executables beside different OpenMM DLLs')
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    cases = [(1000000, 0), (1000000, 750000), (7950256, 5955638)]
    for count, constraints in cases:
        for repeat in range(args.repeats):
            order = ['baseline', 'candidate'] if repeat % 2 == 0 else ['candidate', 'baseline']
            for variant in order:
                exe = binaries[variant]
                p = subprocess.run([str(exe), str(count), str(constraints)], cwd=exe.parent,
                                   capture_output=True, text=True, timeout=240)
                if p.returncode:
                    raise RuntimeError(f'{variant}: {p.returncode}: {p.stderr}')
                row = json.loads(p.stdout)
                row.update(variant=variant, repeat=repeat+1)
                rows.append(row)
                (args.output/'raw.json').write_text(json.dumps(rows, indent=2)+'\n')
                print(json.dumps(row), flush=True)
    summary = []
    for count, constraints in cases:
        selected = [r for r in rows if r['particles'] == count and r['constraints'] == constraints]
        if len({r['input_fingerprint'] for r in selected}) != 1:
            raise RuntimeError('Benchmark inputs differ between DLLs or repeats')
        item = dict(particles=count, constraints=constraints, identical_input_fingerprints=True)
        for variant in ['baseline', 'candidate']:
            a = [r for r in selected if r['variant'] == variant]
            keys = ['construction_seconds', 'destruction_seconds']
            keys += [f'{phase}_{metric}_kib' for phase in ['before', 'constructed', 'destroyed']
                     for metric in ['working_set', 'peak_working_set', 'private']]
            item[variant] = {key: statistics.median(r[key] for r in a) for key in keys}
        summary.append(item)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
