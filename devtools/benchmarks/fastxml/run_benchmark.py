#!/usr/bin/env python3
"""Compare two BenchmarkXmlMemory executables in fresh Linux or Windows processes."""
import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import sys
import statistics
import subprocess

if sys.platform != "win32":
    import resource


def address_space_limit():
    resource.setrlimit(resource.RLIMIT_AS, (6*1024**3, 6*1024**3))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fixtures', type=Path, help='Optional existing stock-generated fixture directory')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    binaries = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
    args.output.mkdir(parents=True, exist_ok=True)
    fixtures = args.fixtures or args.output/'fixtures'
    fixtures.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['OPENMM_FAST_XML'] = '1'
    env['OPENMM_FAST_SYSTEM_WRITE'] = '1'
    env['OPENMM_FAST_STATE_WRITE'] = '1'

    def execute(variant, operation, kind, count, path):
        p = subprocess.run([str(binaries[variant]), operation, kind, str(count), str(path), 'file'],
                           capture_output=True, text=True, timeout=180, env=env,
                           preexec_fn=address_space_limit if sys.platform != 'win32' else None)
        if p.returncode:
            raise RuntimeError(f'{variant} {operation} {kind}: {p.stderr}')
        return json.loads(p.stdout)

    cases = [('system', 500000), ('system-fallback', 500000), ('state-pv', 1000000)]
    expected = {}
    manifest = {}
    for kind, count in cases:
        path = fixtures/f'{kind}-{count}.xml'
        if not path.exists():
            expected[kind] = execute('baseline', 'generate', kind, count, path)['digest']
        with path.open('rb') as f:
            manifest[path.name] = {'bytes': path.stat().st_size,
                                   'sha256': hashlib.file_digest(f, 'sha256').hexdigest()}
    (args.output/'fixtures.json').write_text(json.dumps(manifest, indent=2)+'\n')
    groups = collections.defaultdict(list)
    with (args.output/'results.jsonl').open('w') as results:
        for kind, count in cases:
            for operation in ['read', 'write']:
                for repeat in range(args.repeats):
                    order = ['baseline', 'candidate'] if repeat%2 == 0 else ['candidate', 'baseline']
                    for variant in order:
                        path = fixtures/f'{kind}-{count}.xml' if operation == 'read' else Path(os.devnull)
                        if operation == 'read':
                            with path.open('rb') as f:
                                while f.read(4*1024*1024):
                                    pass
                        row = execute(variant, operation, kind, count, path)
                        expected.setdefault(kind, row['digest'])
                        if row['digest'] != expected[kind]:
                            raise RuntimeError(f'Object digest mismatch: {row}')
                        row.update(variant=variant, operation=operation, kind=kind,
                                   particles=count, repeat=repeat+1)
                        results.write(json.dumps(row)+'\n')
                        results.flush()
                        groups[(kind, operation, variant)].append(row)
                        print(json.dumps(row), flush=True)
    summary = []
    for (kind, operation, variant), rows in groups.items():
        item = dict(kind=kind, operation=operation, variant=variant,
                    particles=rows[0]['particles'], samples=len(rows),
                    peak_scope=rows[0].get('peak_scope', 'operation' if sys.platform != 'win32' else 'process_lifetime'))
        for metric in ['seconds', 'baseline_kib', 'peak_kib', 'retained_kib']:
            values = [r[metric] for r in rows]
            item[metric+'_median'] = statistics.median(values)
            item[metric+'_min'] = min(values)
            item[metric+'_max'] = max(values)
        summary.append(item)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    # Canonical output equality is checked separately from timed operations.
    verification = []
    for kind, count in cases:
        path = args.output/'candidate-output.xml'
        execute('candidate', 'generate', kind, count, path)
        with path.open('rb') as f:
            sha = hashlib.file_digest(f, 'sha256').hexdigest()
        reference = manifest[f'{kind}-{count}.xml']['sha256']
        if sha != reference:
            raise RuntimeError(f'XML byte mismatch: {kind}')
        verification.append({'kind': kind, 'particles': count, 'sha256': sha, 'byte_identical': True})
        path.unlink()
    (args.output/'verification.json').write_text(json.dumps(verification, indent=2)+'\n')


if __name__ == '__main__':
    main()
