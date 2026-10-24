# Four CUDA host-memory and startup/checkpoint optimizations

Base: 0fe8ca2f27f0c38650d9d9ea904d935b00cf9bda, fastxml-experimental-20261024.
All four reviewed changes are separate commits. Exact release commits are in the bundle's validation/release.json. The restored fastxml-clean-v4 branch remains at 655c38e9179d984f095f32d811d3eceed7c69534.

## Implementation

1. VirtualSiteInfo now owns a private immutable compact snapshot: exact RTTI pointers, size_t parent offsets, flat captured indices, and three exact zero-padded double weight slots. Legacy ignored weights, arbitrary parent lists, NaN/Inf/signed-zero identity and snapshot ownership are preserved. No public class layout changes.
2. validateAtomOrder searches already sorted molecule offsets instead of building a node-based set for each relative atom. Existing iteration order, first mismatch, reset and reorder behavior remain.
3. CUDA RNG checkpoint save/load uses at most 8 MiB of host staging per transfer. A private CMake flag selects this for both shared and static CUDA targets. Device byte offsets and total sizes are checked; the owning Context/current stream is selected and synchronized; reads complete before their chunk is uploaded and stream/CUDA errors propagate. randomPos, generated/grown pools, seed bytes, checkpoint format and random-number sequence are unchanged. Other backends retain their original common fallback. Partial failure is not transactional.
4. Atom reordering retains the center/wrapping/bin/sort calculations and records the full permutation. Position, mixed correction and velocity arrays are transferred in separate scopes, with zero padding, original index mapping, cell offsets and final listener notifications preserved. No GPU kernel changes. Partial failure is not transactional.

The core OpenMM.dll is unchanged. All 531 CUDA export names/ordinals and public headers remain unchanged. CUDA static and non-CUDA fallback compilation checks passed; those are compilation evidence, not full static/OpenCL/HIP runtime validation.

## Measurements

RTX 5090, driver 621.01, CUDA 13.4.92, MSVC 19.51/O2/MD, Windows. Each case uses 2,000,000 synthetic particles, one excluded warmup per DLL and three alternating fresh-process trials per DLL. Medians are below. External load varied: earlier correctness tests encountered heavy GPU utilization, while benchmark pre-run snapshots were 0-12% utilization. Per-trial snapshots are retained; timings are scoped observations, not MD/FAH speed claims. Windows process-lifetime working-set peaks are sampled before verification buffers are created. They include earlier construction peaks and must not be interpreted as reset operation-only peaks.

Comparisons isolate cumulative stages: base->01 for snapshots, 01->02 for validation, 02->03 for RNG, 03->final for reordering.

| Change/metric | Baseline | Candidate |
| --- | ---: | ---: |
| Snapshot post-constructor working set, MiB |586.88|545.99|
| Snapshot constructor lifetime peak, MiB |632.66|592.29|
| Snapshot Context construction, seconds |2.4549|2.1527|
| Five atom-order validations, seconds |1.15946|0.30839|
| RNG save lifetime peak, MiB |711.04|633.79|
| RNG load lifetime peak, MiB |711.36|633.81|
| Reorder scene lifetime peak, MiB |1251.16|1250.72|

The measured snapshot retained-working-set reduction is40.89MiB at this size. This is post-constructor host memory, not a measured long-running MD average or VRAM saving. Validation is faster, but its lifetime peak is masked by construction (592.8 versus592.6MiB). RNG save/load peak reductions are about 77 MiB in these bounded I/O tests; the candidate's earlier constructor highwater masks some of the staging reduction. Do not add savings from different phases.

Startup/validation scenes contain particles, constraints and virtual sites, with no Force objects and no integration. RNG scenes have no force and explicitly initialize/generate the default RNG pool. Save writes to Windows NUL; load supplies a synthetic bounded zero stream and verifies byte consumption. Their measured save/load times are0.01515->0.00864s and0.01726->0.00934s. These are isolated RNG serialization timings, not disk, compression, full-checkpoint or FAH wrap-up timings. Exact generated-state correctness is tested separately.

The reorder scene uses a periodic cutoff NonbondedForce with zero epsilon/charge and excludes within-water pairs. Its measured operation is setPositions+computeVirtualSites+three forced reorders. Times0.39099->0.38124s have overlapping ranges; no meaningful reorder speedup is established. Lifetime peak is essentially unchanged because initialization already peaked higher.

### Reorder phase observation

A separate caller with start/end stdin acknowledgments lets an external PSAPI sampler measure only that active phase on the same uninstrumented DLLs. The caller blocks after the end marker until sampling is stopped and joined, then verification may allocate memory. Three alternating trials per DLL, with excluded warmups, produced median sampled phase peaks1227.625->1098.926MiB, about128.70MiB lower. Baseline range 1218.84-1232.54MiB; candidate 1098.16-1110.32MiB. These sampled maxima are observed lower bounds on actual phase peaks; their difference is not a lower bound on exact savings. Requested1ms sampling yielded actual median intervals near 15.6 ms on Windows, with21-27samples per measured trial. Raw timestamps, working sets, private bytes, handshakes and caller hashes are retained in reorder-phase/. Original lifetime records remain in benchmarks/.

## Correctness and failure history

* Snapshot actual-helper oracle:2 final host passes, also independently rerun. Covers known/derived/custom RTTI, arbitrary 1000-parent lists, local-coordinate ignored weights, nonfinite/signed-zero values and immutable lifetime across System replacement/destruction.
* Final regression records:54passes, comprising24selected existing CUDA test successes,12controlled integrator suite successes and18molecule/reorder/checkpoint cases. The latter include all three precisions, both checkpoint directions and deliberately invalid index maps repaired without changing logical positions/velocities. Existing triclinic PME reordering assertions are retained.
* One original unseeded BASELINE LangevinMiddle single test failed its explicitly stochastic initial-temperature assertion:303.873K versus300K. It is preserved in regressions/raw.json and was not retried to obtain a pass. The controlled scratch suite keeps all assertions/steps/tolerances, adds default integrator seeds20261024 and fixes the initial-temperature velocity seed; explicit5/10 seed tests still override it. All six controlled baseline/candidate cases pass. One separate massless-particle velocity initialization remains unseeded, so the entire scratch suite is not claimed fully deterministic. Controlled source/diff is supplied.
* RNG 15/15: generated/grown pools beyond two chunks, byte-identical saved full checkpoints in all precisions (direct comparison and SHA256 in rng/exact-bytes.json), self-replay and cross-DLL continuation both directions. Actual private helper tests cover1,8MiB+/-1,and two-chunk+tail sizes; short/exception-enabled reads/writes; incomplete chunks not uploaded; non-default queued transfers; and genuinely uninitialized RNG. A helper stream synchronization test does not assert unrelated full-state fields become coherent on arbitrary custom queues.
* Nonperiodic cutoff reordering:6/6baseline/candidate cases,260integration steps and exact finite state/atom-index hashes and self-checkpoint continuation in all precisions.
* Fresh-cache isolated PME succeeds with CUDA environment variables removed and OS-only PATH. All five OpenMM/plugin/NVIDIA modules load from the staged package. Core and export checks are included.

Initial stage-smoke.json is historical incremental evidence before the parent-count robustness cast; final combined tests and final helper reruns validate the delivered implementation. Independent final review and source/binary identities accompany this report.

## Limits

No live FahCore work unit or installation was changed or executed. No PPD, actual work-unit end-to-end startup/wrap-up, VRAM saving, reduced precision or modified RNG sequence is claimed. XML formatting/compression are unchanged. The pre-existing ComputeForceInfo destruction/ownership issue is deferred; no destructor/vtable ABI change or leak repair is included. Measurements concern the first live Context and selected operations. No exact additive RAM estimate should be extrapolated to the7.95million-particle FAH topology.
