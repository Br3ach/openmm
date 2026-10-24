# Context State collection on Windows

Baseline: `85570f1511100166effb157eba6f5f72ff7c9693`, fastxml-clean-v4.
Candidate: ownership transfers within `Context::getState()` only. The library
uses its existing StateBuilder rvalue setters and `takeState()`. Its original
copying setters and `getState()` remain unchanged for other callers.

## Why this path matters

Static analysis of the supplied Windows FahCore_28.exe (SHA-256
`8f8ff4db1916f8cb3400edd838b69b5befaceea55518ae5de83d17629b5bb815`)
finds `getState(0x1f, false, -1)` before periodic checkpoint validation
(call at image VA `0x140029a2a`) and during normal unencrypted completion
(`0x14002a870`). The former is followed by binary `createCheckpoint()`
(`0x14002a010`); the latter by State XML serialization (`0x14002a923`).
The executable is not included here. This is static inspection, not a run of
the modified library inside the FAH core.

Unlike the earlier XML-only changes, this change reaches State collection at
regular checkpoint boundaries. It changes neither binary checkpoint bytes nor
GPU download/kernel work. It requires replacing the core OpenMM DLL; no caller
recompile or platform-plugin rebuild is needed for this internal change.

## Measurements

Three fresh processes per variant/case, alternating order, one identical caller
executable with different adjacent OpenMM DLLs. Median values follow. Mask 3
requests positions and velocities. Mask 31 matches the core's full State request
(positions, velocities, forces, energy, parameters). There are no forces or
constraints in the synthetic System, so force arrays are zero and kinetic
energy is computed from deterministic velocities.

| Particles | Mask | Before seconds | After seconds | Before peak MiB | After peak MiB |
|---|---|---|---|---|---|
| 1,000,000 | 3 | 0.017574 | 0.006382 | 204.0 | 158.3 |
| 1,000,000 | 31 | 0.032237 | 0.017914 | 249.9 | 181.2 |
| 7,950,256 | 3 | 0.167743 | 0.056981 | 1584.3 | 1220.3 |
| 7,950,256 | 31 | 0.286854 | 0.146588 | 1948.2 | 1402.3 |

All returned-array and energy FNV-1a fingerprints match across every run. The
benchmark measures one `getState()` call, excluding construction and hashing.
Peak RAM is the process-lifetime working-set high-water mark, including setup.
The full-State 7,950,256-particle case drops about 546 MiB of transient peak,
consistent with removing a duplicate of three double-precision Vec3 arrays.
Retained State payload memory is essentially unchanged. Shared-desktop timing
variation remains possible.

These are Reference-platform State-collection measurements. They do not
measure CUDA transfers, real FAH force evaluation, binary disk-write time,
XML/BZip2 output, end-to-end checkpoint pauses, integration throughput or PPD.
The Reference baseline memory and timings are not predictions of CUDA totals.

## Correctness

43 tests pass: all 39 serialization executables, the existing Reference
checkpoint, CustomIntegrator and kinetic-energy suites, and the new State
snapshot regression. The latter exercises all 128 data-flag combinations,
force-group filtering, periodic wrapping without Context mutation, and retained
State values/XML after later retrievals, Context destruction, and changes to
positions, velocities, box, parameters, time, step count and integrator data.

`Context::getState()` still performs the same calculations and collects the
same data in the same order. Only ownership of method-local buffers is
transferred. Returned States continue owning independent arrays and metadata.

## Reproduce

Build the baseline and candidate core DLLs with the same x64 Release options.
Compile `BenchmarkContextState.cpp` once with the baseline headers/library:

```powershell
cl /nologo /O2 /Ob2 /MD /EHsc /std:c++14 /DWIN32 /D_WINDOWS /I"$sourceDir/openmmapi/include" /I"$sourceDir/olla/include" /I"$sourceDir/serialization/include" /Fe:BenchmarkContextState.exe "$sourceDir/devtools/benchmarks/fastxml/BenchmarkContextState.cpp" "$baselineBuild/OpenMM.lib" Psapi.lib
```

Copy that identical executable beside each DLL, then run:

```powershell
python "$sourceDir/devtools/benchmarks/fastxml/run_context_state.py" --baseline "$baselineBuild/BenchmarkContextState.exe" --candidate "$candidateBuild/BenchmarkContextState.exe" --output "$resultsDir"
```

The Windows harness requires a 6 GiB per-process committed-memory limit via a
Job Object and aborts if it cannot establish it. The runner applies a 120-second
timeout per process. Raw runs, medians, compiler/source identity and test
outputs are retained in this directory.
