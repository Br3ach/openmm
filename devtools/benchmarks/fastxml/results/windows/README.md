# Native Windows validation

Validated runtime commit 617b8051 against parent 655c38e9 on 24 October 2026.
The follow-up fixes a Windows regression-test linker failure by delegating to the registered State proxy instead of subclassing the unexported StateProxy. Runtime library sources are unchanged.

39/39 serialization tests passed. All three checked synthetic outputs are byte-identical between versions.

## Synthetic workloads

Medians of three fresh-process runs. Peak is the Windows lifetime working-set high-water mark in MiB, including setup allocations. It is not directly comparable to Linux operation-reset RSS.

| Fixture | Operation | Baseline peak MiB | Candidate peak MiB | Baseline seconds | Candidate seconds |
|---|---|---:|---:|---:|---:|
| system | read | 79.5 | 79.6 | 0.534 | 0.517 |
| system | write | 1379.5 | 77.9 | 8.644 | 0.368 |
| system-fallback | read | 1706.3 | 79.6 | 11.544 | 0.446 |
| system-fallback | write | 1679.4 | 77.8 | 8.746 | 0.320 |
| state-pv | read | 54.4 | 54.6 | 0.788 | 0.734 |
| state-pv | write | 119.8 | 119.8 | 0.911 | 0.911 |

System fixtures have 500,000 particles. The State control has positions and velocities for one million particles. Raw runs, fixture fingerprints, ranges and output comparisons are in synthetic/.

## Actual Folding@home inputs

Used stable copies from project 18289, run 319, clone 1, generation 0: 7,950,256 particles, 5,955,638 constraints and six forces. System XML is 1,927,759,124 bytes; State XML is 1,947,913,030 bytes and includes positions, velocities, forces and energies. State step count: 599,989. XML and binary checkpoints are not committed.

Medians of three fresh-process reads, with alternating variant order and warm files:

| Input | Baseline peak MiB | Candidate peak MiB | Baseline seconds | Candidate seconds |
|---|---:|---:|---:|---:|
| system | 1500.2 | 1500.0 | 15.579 | 15.966 |
| state | 554.4 | 554.4 | 21.332 | 21.134 |

This System already fits the old fast reader. The patch does not materially reduce its loading memory. Timing variation is visible in raw runs; these samples do not establish a sustained throughput gain.

Both complete candidate round trips match their original XML byte-for-byte after changing only the root openmmVersion attribute from 8.2 to 8.3.1. SHA-256 fingerprints and separate read/write measurements are in fah/verification.json. The two-byte difference is expected when serializing older inputs with the newer library.

The actual System was not written through the original whole-System tree. Synthetic comparisons establish that result while avoiding an unnecessarily large legacy allocation on a working FAH machine. The real State comparison measures reads; its candidate output is verified separately.

## Existing callers

An unchanged executable compiled against the old header and library loads successfully with the candidate DLL. Its late-fallback read benefits from the new reader, while its System writer retains the original tree path. See old-header-caller.json. Recompile callers to obtain the new System writer.

## Scope and reproduction

See the parent [benchmark README](../../README.md) for Linux and Windows commands. Compiler flags, source hashes and the direct native build method are recorded in source-and-build.json.

These runs validate native serialization and real work-unit XML compatibility. CUDA Context initialization, a rebuilt FAH core, trajectory correctness, PPD, binary checkpoints, BZip2 cost and OS memory compression remain outside this validation. The Linux benchmark changes were not re-executed in this Windows session.
