# Earlier release of constraint-validation storage

Baseline: `ab4179b71b83e876bd6b5c0888137da268bc3d01`. The only runtime change
adds an inner scope around the existing `constraintAtoms` set and validation
loop in the `ContextImpl` constructor. The set is unused after that loop. Its
allocation, validation order, `find` plus `insert` calls, error messages and
public ABI remain unchanged.

Previously the tree survived property validation, ForceImpl creation and
`Platform::contextCreated()` / `linkedContextCreated()`. Those hooks construct
Reference constraint data or CUDA IntegrationUtilities, including large host
constraint/SETTLE/SHAKE/CCMA preprocessing structures. The new scope releases
the validation tree before those allocations. The old code already destroyed
it before `ContextImpl::initialize()`: this change does not remove overlap
with every later kernel initialization, nor change integration/checkpoint work.

## Constructor A/B measurements

Three fresh processes per variant and case, alternating order, using a byte-
identical Windows caller beside the two DLLs. The synthetic System has positive
particle masses, mostly disjoint three-constraint triangles with water-like
distances, at most one incomplete triangle, and unconstrained remaining atoms.
It has no forces or virtual sites and does not represent a physical simulation.
The largest case matches the particle and constraint counts scanned from the
available project 18289 / run 319 / clone 1 / generation 0 System XML, but does
not load or construct that actual FAH System.

Timing covers the full public Reference Context constructor, including
`ContextImpl::initialize()`. System creation, input hashing and Context
destruction are excluded from construction time. Destruction is timed separately.
The no-constraint case controls for unrelated changes. Every run verifies the
constructed Context and matching input fingerprints.

| Particles | Constraints | Before constructor seconds | After constructor seconds | Before lifetime peak MiB | After lifetime peak MiB |
|---|---|---|---|---|---|
| 1,000,000 | 0 | 0.079151 | 0.080046 | 154.5 | 154.6 |
| 1,000,000 | 750,000 | 0.356466 | 0.346357 | 304.9 | 270.7 |
| 7,950,256 | 5,955,638 | 3.331123 | 3.245343 | 2987.2 | 2713.7 |

At 7,950,256 particles and 5,955,638 constraints, the observed process-lifetime
peak falls by 273.4 MiB (about 9.2%). Unlike the earlier isolated set-allocation
probe, this is a measured full Reference Context construction result. It is
still synthetic and does not establish the CUDA/FAH process peak reduction.
Current working set after construction is about 1,066.9 vs 1,063.8 MiB, so the
benefit is transient overlap, not a sustained 273 MiB reduction. Do not add this
peak difference to different-phase getState/getMolecules savings to predict a
whole-work-unit peak.

Raw records separately retain current working set, lifetime peak working set
and current private committed memory before construction, after construction
and after destruction. No private-memory peak is measured. Windows lifetime
peaks include System setup and cannot be reset; here the construction phase
sets the high-water mark. No EmptyWorkingSet/heap trim is performed. Timing
varies on a shared desktop; the small time difference does not establish a
general CPU-throughput, checkpoint or PPD improvement.

## Correctness and stochastic-test triage

All 50 candidate test executables were run: the previous 49 plus the new
`TestContextConstraintValidation`. Initially 49 passed; the existing Reference
CustomIntegrator suite failed its explicitly stochastic initial-temperature
assertion (296.618 K versus 300 K, with 1% tolerance). This is preserved in
`initial-tests.json`. Three further fresh baseline/candidate runs are recorded
in `stochastic-retries.json`: baseline passed 3/3, candidate 1/3. These are not
reported as an unconditional 50/50 clean run.

A scratch build of the full CustomIntegrator suite changed only the affected
`setVelocitiesToTemperature(targetTemperature)` call to supply seed 20261024.
The same scratch executable passed with both DLLs (`fixed-seed-test.json`).
Repository test sources and tolerances were not changed.

To control this independently of the chosen seed, `CompareSeededTemperature.cpp`
constructs a zero-constraint 50,000-particle Reference Context and compares seeds
0..99 plus 20261024. All 101 temperatures and fingerprints of every velocity
component's IEEE bytes match exactly across the two DLLs. Eight seeds produce
out-of-1% temperatures with BOTH DLLs (24, 33, 40, 50, 52, 53, 69, 83), showing
that this temperature tolerance can fail independently of the scope change.
The full raw comparison and summary are retained. The controlled full suite
and 101-seed comparison resolve the observed stochastic failures without
weakening the test or modifying production numerical code.

The new regression passes unchanged against both DLLs. It covers no constraints,
valid triangles, constraints between massless particles, Context reinitialization,
self-pairs, negative/high indices, mixed massless/massive constraints, duplicate
and reversed pairs, and precedence between competing validation errors and
invalid platform properties. It checks exact existing exception strings.

## Independent review and reproduction

A separate reviewer approved the scope-only design before implementation.
Final review evidence is recorded in `review.json` after approval.

Build both core DLLs with matching x64 Release options. Compile
`BenchmarkContextConstruction.cpp` once against the matching import library
and headers, link Psapi.lib, and place identical copies beside both DLLs:

```powershell
python "$benchmarkSource/run_context_construction.py" --baseline "$baselineBuild/BenchmarkContextConstruction.exe" --candidate "$candidateBuild/BenchmarkContextConstruction.exe" --output "$resultsDir"
```

The benchmark requires a 6 GiB per-process committed-memory Job Object limit;
the runner enforces a 240-second timeout. To reproduce stochastic triage, build
`CompareSeededTemperature.cpp` once with the same headers/import library and
run it beside both DLLs, comparing every seed/temperature/fingerprint record.
For the scratch full suite, copy `TestCustomIntegrator.h` into a scratch include
directory, supply seed 20261024 only at the initial-temperature call, and compile
`TestReferenceCustomIntegrator.cpp` using that scratch header. Preserve the
existing tolerance and all other test code. Executable and DLL hashes are in
`build.json`.
