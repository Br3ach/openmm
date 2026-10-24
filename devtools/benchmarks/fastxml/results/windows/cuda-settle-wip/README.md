# CUDA SETTLE work in progress, preserved during branch separation

This is an unfinished snapshot, not a released CUDA plugin. At the user's
request, fastxml-clean-v4 is being restored to its pre-24 October revision,
655c38e9179d984f095f32d811d3eceed7c69534. All five commits made on 24 October
through 19c091248e636747d74d516e28c795e73d9063a2 are preserved on
fastxml-experimental-20261024, together with this later uncommitted work.

## Code retained

- CUDA 13 context-creation compatibility: new four-argument call with NULL
  parameters for CUDA_VERSION>=13000, original call for older SDKs. This
  compatibility change is separated from the SETTLE optimization commit.
- Private 20-byte, sorted, at-most-two-neighbor SETTLE records in place of
  temporary maps. Preserve float conversions, clear cascade, solver-selection
  branches, int4(w=0)/float2 upload ordering and existing public API/layout.
- Release temporary records/clusters after SETTLE uploads, before later
  SHAKE/CCMA/virtual-site setup. Other solver and kernel code is unchanged.
- Host equivalence regression using the exact production helper, full-topology
  host benchmark, fixed-seed mixed-solver CUDA snapshot/checkpoint test,
  bounded synthetic CUDA startup benchmark and isolated runtime-loader probe.

## Validation actually completed

The same core DLL from 19c0912 was used in both builds. The baseline plugin
contains only the shared CUDA 13 compatibility fix; the candidate reuses all
baseline objects except IntegrationUtilities.cpp. No NVIDIA GPU math kernels
were modified.

- TestSettleConnectivity passed against both build directories. A separate
  reviewer independently reran it and both full-topology host modes.
- Six fresh host processes using actual 18289 masses/constraints matched every
  ordered SETTLE output and remaining-neighbor fingerprint. They identify
  1,984,149 clusters from 7,950,256 particles and 5,955,638 constraints. These are
  host-model measurements, not CUDA Context/FAH peak measurements. Raw results
  and all medians are in host/. Instrumentation/model limitations are the same
  as the earlier assessment: output tuples 20B versus production int4+float2
  24B, entry-count audit/memory queries in timing, and warmed row destruction.
- Baseline CUDA SETTLE single-precision smoke passed.
- Candidate completed 16 CUDA suite/precision cases: 15 passed. The original
  CustomIntegrator single-precision test failed its documented stochastic
  initial-temperature tolerance (296.361 K versus 300 K). The mixed-precision
  CustomIntegrator case passed. Original outcomes are retained unchanged.
- The next candidate NoseHoover mixed-precision suite timed out after 600 seconds.
  The runner stopped there. The remaining mixed Monte Carlo and all double-
  precision cases were not run. The timeout has NOT been attributed to the
  optimization or dismissed as pre-existing: a baseline comparison is pending.
- A fixed-seed scratch CustomIntegrator caller was built but not run. Its
  initial-temperature call alone supplies seed 20261024; test tolerances and
  repository numerical test headers remain unchanged.

## Validation still outstanding

Baseline comparison/diagnosis of the NoseHoover timeout, controlled-seed followup
for the stochastic failure, remaining CUDA precision/suite cases, matched-seed
baseline/candidate snapshots, cross-DLL binary checkpoint continuation, scaled
CUDA startup memory measurement, and isolated runtime-bundle loading/PME/FFT
validation are incomplete. Snapshot, startup and loader utilities were compiled
but had not been executed when work was redirected. No CUDA memory or PPD gain
has been established for this implementation, and no final plugin package was
published.

The reviewer approved the design before edits and found no early code blocker;
final approval here is only for preserving this clearly labelled WIP snapshot.
It is not production-release readiness approval. See review.json once recorded.

## Dependencies and reproduction notes

The CUDA 13.4 build imports cufft64_12.dll and nvrtc64_130_0.dll; NVRTC also needs
nvrtc-builtins64_134.dll from toolkit bin/x64. These three files are not included
in this source snapshot. An older FAH runtime bundle cannot be assumed to supply
them. A verified bundle with matching runtimes/licenses was planned but not yet
prepared. Do not treat this as a self-contained plugin replacement.

workspace-scripts/ preserves the original local build/test drivers. Their paths
and expected workspace layout are recorded as WIP provenance, not a portable
installer. build-config.json retains exact compiler flags, source list and link
inputs. No DLLs, CUDA runtimes, original FAH XML or derived 159 MB constraint binary
are committed. Local artifacts remain in the workspace for continued work on
the experimental branch.
