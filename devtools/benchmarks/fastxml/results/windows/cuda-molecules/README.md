# CUDA molecule-connectivity memory optimization

Baseline: `d9660c6c0c8edfd42ff62e0125eee38e53eb0c21` on `fastxml-experimental-20261024`, including earlier CUDA 13 compatibility and compact SETTLE preparation. Stable `fastxml-clean-v4` remains at `655c38e9179d984f095f32d811d3eceed7c69534`.

## Change

`ComputeContext::findMoleculeGroups()` joins constraint and ComputeForceInfo group-chain connectivity directly using temporary union/find storage. Components retain minimum-atom order and ascending members. DSU scratch is destroyed before descriptors are populated. Atom rows and exhausted instance/offset buffers move into persistent owners; representative indices replace deep copies of Molecule objects. All constraints remain connectivity edges, including massless ones. ForceInfo-only connections cannot be replaced with core getMolecules connectivity.

VirtualSiteInfo registration, both group-enumeration callback passes, comparison loops, first representative selection, invalidation cache behavior, public headers/layout, kernels, solver selection and checkpoint format retain their prior behavior. Only ComputeContext.cpp and its new noninstalled MoleculeConnectivity.h change runtime code. Core OpenMM.dll is unchanged. This common host code is validated here through CUDA; OpenCL/HIP are untested.

## Measured memory

Windows PeakWorkingSetSize is a process-lifetime peak sampled after the complete CUDA Context constructor. These are actual uninstrumented baseline/candidate CUDA DLLs, serial fresh processes, mixed precision. These startup scenes contain no Force objects and perform no integration; water scenes add constraints and virtual sites. Each case has one excluded warmup per DLL followed by three alternating measured trials per DLL. Values below are medians. Drivers/toolkit/hardware match across variants: RTX 5090, driver 621.01, CUDA 13.4.92, MSVC 19.51, /O2 /MD.

| Synthetic scene | Baseline peak MiB | Candidate peak MiB | Reduction MiB |
| --- | ---: | ---: | ---: |
| 250,000 unconstrained particles | 295.91 | 242.55 | 53.36 |
| 250,000 four-site water particles | 310.32 | 241.41 | 68.91 |
| 500,000 four-site water particles | 441.30 | 291.57 | 149.72 |

At 500k the measured peak reduction is 33.93%. Retained working set after construction is effectively unchanged: 280.63 versus 279.92 MiB. PrivateUsage is separately recorded and must not be confused with peak working set. Constructor wall-time medians were 1.5088 versus 1.2489 seconds at 500k; these are descriptive setup timings, not an MD throughput/PPD result. No VRAM saving has been measured.

A separate instrumented baseline scratch build brackets the old highwater within findMoleculeGroups graph/metadata overlap, including the interval between atom copying and staging release. It is not a uniquely attributed single-allocation trace. Instrumented binaries were excluded from performance measurements and delivery.

### Large synthetic outcome

One bounded trial per variant used 7,950,256 particles and 5,962,692 constraints in repeated four-site waters. This matches the FAH fixture's particle count only; its topology and constraint count differ. Each launch required >=8 GiB free VRAM and >=12 GiB available physical host RAM, with a 6 GiB Windows Job process-memory cap and a 180-second child timeout. Both resource gates passed. The baseline failed with `bad allocation` while running under that cap; the error does not prove the cap was the exact cause or establish a measured baseline peak. The candidate completed: process peak working set 1,995,296 KiB (1.90 GiB), retained working set 1,808,960 KiB, PrivateUsage 3,947,016 KiB, Context construction 5.9523 seconds. This is a single synthetic result with no exact paired savings or speedup estimate. No caps were raised and no baseline retry substituted for the failure. Owned children exited; other processes were untouched.

## Correctness validation

* The actual private DSU helper matches the untouched legacy adjacency/findMolecules oracle for empty/singleton/self/repeated/reversed edges, group chains, 1,000 randomized orientations, a 200k chain, nested virtual-site connectivity, repeated enumeration and snapshot ownership. Both DLL directories pass; the reviewer also reran it independently.
* Final finite-checked reordering/equivalence suite: 12 successful fresh processes across single/mixed/double, each DLL saving and self-replaying checkpoints and each loading the other DLL's checkpoint. Full position/velocity/force/energy/time/step hashes, atom-index permutations, and globally interleaved group-enumeration callback traces match exactly. Each scene forces nonidentity periodic cutoff reordering, integrates 260 steps, exercises ForceInfo-only cross-molecule connectivity, two/three/nested virtual sites, four explicit particle/group invalidations, merges, and nonbonded updateParametersInContext. Regrouping retains ForceInfo cache size. Nonfinite states are rejected before hashing. Parallel identity callbacks are not traced; their comparison loops were checked in review.
* 24 original/selected CUDA test processes pass: original Nonbonded reordering/parameter changes/offsets/periodic/exceptions assertions, full VirtualSites, full Checkpoints, and mixed SETTLE snapshots, both DLLs and all three precisions. Existing tolerances are unchanged; SETTLE snapshot fingerprints match across DLLs.
* The first new fixture produced NaNs on the unchanged baseline. Fixture fixes used equal peripheral masses within SETTLE triangles and added nested-site parent exclusions before final reruns. This was not treated as a production regression or hidden by relaxed tolerances.
* All 531 exported CUDA function names and ordinals match. Identical caller executables are used per comparison; hashes are in binaries.json.
* Candidate fresh-cache PME/JIT succeeds with CUDA environment variables removed and only Windows directories on PATH. Module checks confirm core/plugin and all three NVIDIA runtime DLLs load from the staged release directory.

Raw trials, resource gates, failures and binary hashes accompany this report. No live FahCore_28 execution, actual work-unit peak RAM, steady MD performance or PPD is claimed. The package retains prior validation under validation/ and adds this evidence under validation/cuda-molecules/.
