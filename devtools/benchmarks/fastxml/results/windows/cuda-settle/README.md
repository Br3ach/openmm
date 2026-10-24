# CUDA SETTLE storage: final Windows validation

Implementation is on `fastxml-experimental-20261024`. The separate original
`fastxml-clean-v4` branch remains restored to 655c38e. Core OpenMM.dll is unchanged
from 19c0912 and retains all earlier XML/State/molecule/validation-lifetime fixes.
The CUDA 13 context-creation compatibility fix is isolated in 239c69b. The SETTLE
implementation was preserved as WIP in ca68de7; this record completes validation
of the same runtime code and binary rather than changing GPU numerical kernels.

## Runtime change and scope

A private noninstalled 20-byte sorted two-neighbor record replaces each temporary
std::map used by IntegrationUtilities SETTLE discovery. The active-degree 2
filter bounds eligible neighbors to two. Preserve duplicate-key updates,
double-to-float conversion before comparisons, ascending neighbor order, the
ascending clear cascade, central-atom branch order, isShakeAtom flags and exact
int4(w=0)/float2 upload order. The records and temporary cluster list now leave
scope immediately after SETTLE uploads, before SHAKE/CCMA/virtual-site setup.
Public headers/layout, device arrays, solver formulas, kernels and checkpoint
format are unchanged. Existing caller/plugin API exports match: 531 names and
ordinals. The core DLL itself is byte-identical in baseline and candidate.

CUDA, OpenCL and HIP compile this common source into separate platform libraries.
Only Windows CUDA is built/validated here. No OpenCL/HIP binary is supplied.

## Actual CUDA Context memory benchmark

Same caller executable, same core, and matching CUDA 13.4 baseline/candidate
plugins. Only IntegrationUtilities.cpp is recompiled for the candidate; 32 other
plugin translation units are reused. Three fresh processes per variant/case,
alternating order, following one excluded warmup for each variant/case. Mixed
precision, 6 GiB process committed-memory limit and 180-second timeout. Synthetic
water-like cases have one three-constraint triangle plus a massless three-
particle-average virtual site per 4 particles and no forces. This is not an
actual FAH System or MD trajectory benchmark.

The caller subclasses CudaPlatform to sample memory immediately after the
contextCreated hook returns, before ContextImpl::initialize, and samples again
after full public Context construction. Windows PeakWorkingSetSize is a lifetime
high-water mark, not an operation-only allocation total; before/after current
working set and current private commit are retained separately in raw.json.

| Particles | Case | Before through platform creation MiB | After through platform creation MiB | Before whole Context peak MiB | After whole Context peak MiB |
|---|---|---|---|---|---|
| 250,000 | unconstrained | 229.7 | 214.3 | 296.1 | 296.0 |
| 250,000 | water | 263.4 | 230.8 | 310.3 | 310.2 |
| 500,000 | water | 322.9 | 256.9 | 441.7 | 441.3 |

At 500,000 particles, lifetime peak through CUDA platform creation falls by about
66 MiB. The final whole-Context peak remains essentially unchanged, because a
later initialization stage sets a higher high-water mark. This benchmark does
not isolate which later allocation sets it. Do not report 66 MiB as a reduction
in the overall process peak or in retained/VRAM usage. Measured median constructor
times are in summary.json, but these small synthetic timings do not establish
FAH startup/MD/checkpoint/PPD improvements.

Recorded GPU status during these measured runs was around 30 GiB free and 0-8%
utilization. Earlier correctness tests ran under varying external load. The
benchmark is explicitly bounded to 500,000 particles. No full-size CUDA Context
or actual FAH core was run. Full-size host evidence below must not be called a
CUDA/FAH process-memory result or added to different-phase memory savings.

## Full-topology host equivalence

BenchmarkSettleStorage uses the ACTUAL private production helper against std::map.
The extracted real 18289 masses/constraints contain 7,950,256 particles,
5,955,638 active constraints, 11,904,894 eligible neighbor entries and 1,984,149
SETTLE clusters. All six fresh-process ordered-output and surviving-neighbor
fingerprints match. Requested container storage is 878.7 -> 151.6 MiB; model process
peak is 1383.1 -> 498.0 MiB. This is a host preprocessing model, not IntegrationUtilities
or a complete CUDA Context. It omits device work and later setup; its output tuple
is 20B versus production's int4+float2=24B, its shared flags allocate slightly
earlier, timings include audit/memory queries, and row destruction follows an
untimed fingerprint pass that warms the rows. Full raw/median evidence is in host/.

TestSettleConnectivity exercises the actual helper, comparing complete ordered
records/parameters/flags with the legacy map algorithm on 500 deterministic random
graphs, 108 float-sensitive triangles, a 20,000-particle degree 2 cycle, open paths,
4/5/12-cycles, reversed edges, duplicate keys, massless cases and overflow/const
lookup behavior. It passed, including independent reviewer reruns.

## GPU correctness, stochastic test and timeout

Nine existing candidate suites were exercised in single/mixed/double precision:
SETTLE, virtual sites, checkpoints, Verlet, LangevinMiddle, CustomIntegrator,
VariableVerlet, NoseHoover and MonteCarloBarostat. The matrix has 26 direct/recheck
passes and one original stochastic temperature failure. Do not describe this as
an unconditional 27/27 original run.

The original CustomIntegrator single-precision case reported 296.361 K versus 300 K
with its documented 1% stochastic tolerance. A scratch copy of the complete suite
supplied only fixed seed 20261024 to the initial-temperature call; the same
scratch executable passed both baseline and candidate. Repository numerical test
headers, step counts and tolerances were not changed. Original results remain.

An initial candidate NoseHoover mixed case timed out at 600 seconds. It was then
rerun unchanged against BOTH plugins with the same 600-second cap and cache/PATH:
baseline passed 103.53s, candidate 58.59s. The original timeout is not reproduced;
its cause remains undetermined. Do not claim workload caused it or use these
varying-load times as a speedup. Owned-child CPU/memory samples and all original
outcomes are retained. All 26 subsequently added baseline/candidate cases pass,
including baseline SETTLE/vsites/checkpoints, all candidate double cases, the
controlled CustomIntegrator pair and baseline/candidate FFT in all precisions.
FFT results are checked against the existing CPU reference/tolerances.

TestCudaSettleSnapshots is a small fixed-seed scene with all central-atom ID
positions, reversed edges/float-equal lengths, SETTLE, a SHAKE star, scalene/open-
chain/cyclic CCMA constraints, virtual-site forces and full State fingerprints.
Across 12 successful processes, all three precision fingerprints match between
plugins. Save cases verify internal checkpoint replay; each plugin also loads
the other's binary checkpoint and reproduces continuation in both directions.
The runner owns cross-plugin resume verification; load mode correctly does not
claim its own internal replay check. This validates these two CUDA 13.4 builds,
not arbitrary older SDK binaries or execution inside FahCore_28.exe.

## Isolated runtime dependencies and delivery

The CUDA 13.4 plugin needs cufft64_12.dll, nvrtc64_130_0.dll and
nvrtc-builtins64_134.dll. Both isolated baseline/candidate directories passed a
fresh-JIT tiny PME energy check with PATH only Windows/System32+Windows and CUDA
path environment variables removed. GetModuleFileName verified the core, plugin
and all three runtime DLLs loaded from the isolated directory, and recorded hashes
match the supplied files. CUDA toolkit installation/PATH is therefore not needed
for this tested bundle. The installed NVIDIA driver and Microsoft Visual C++ x64
runtime remain prerequisites; tested GPU/driver are RTX 5090/621.01. Driver and
older-GPU compatibility beyond this tested setup has not been established.

The release bundle includes the paired core/plugin, matching import libraries
and public headers, these three unchanged NVIDIA runtime DLLs, OpenMM/third-party
notices and NVIDIA toolkit EULA/LICENSE. Put OpenMMCUDA.dll where the application
loads its CUDA plugin, and place OpenMM.dll plus the three runtime DLLs beside
the application executable. Do this with the application closed and retain its
previous files. Do not change a running FAH core's files. No live application
installation was performed by this validation.

The prior streaming System XML writer requires caller recompilation against its
updated header; the internal State/molecule/validation/CUDA changes are reached
through their compatible DLLs. Actual execution in FAH, FAH peak RAM and PPD
remain unmeasured.

## Independent review

A separate reviewer approved the design before edits, including CUDA 13 compatibility
and the private helper. They independently reran the host regression/full-topology
modes and selected GPU snapshots/SETTLE/checkpoints/isolated PME, verified raw
medians, caller/DLL hashes, exports, dependency paths and reporting limits. Final
release approval is recorded separately in review.json after package review.

All raw outcomes, medians, warmups, original failure and timeout evidence, fixed-
seed results, cross-plugin comparisons and local reproduction scripts are kept
here. The earlier cuda-settle-wip directory remains an historical snapshot of
the intentionally unfinished state during branch reorganization.
