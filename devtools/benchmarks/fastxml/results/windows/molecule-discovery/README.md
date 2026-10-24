# Molecule discovery: native Windows validation

Baseline: `73228c31055cc6fae3ed8f94a0ee33e4f8444c44`, including the earlier
Context State ownership transfer. The only additional runtime change is in
`ContextImpl::getMolecules()`. It uses a method-local disjoint-set forest to
join constraint, force and virtual-site edges instead of retaining a combined
bond vector, per-particle adjacency vectors and traversal buffers.

The exported static `findMolecules()` algorithm, public headers/layout and
callback ordering remain unchanged. The result still orders molecules by their
smallest particle and orders members ascending. It is cached only after complete
construction. Existing callers and platform plugins use the revised core DLL
without recompilation. The previous State-collection improvement is retained.

## Measurements

Three fresh processes per variant/case, alternating order, one byte-identical
caller executable beside two different DLLs. The Reference Context is created
before timing. A custom no-kernel ForceImpl supplies synthetic connectivity:
water-like triangles, a single chain, or isolated particles. The benchmark does
not simulate physical water or load actual FAH input. Sizes include the same
7,950,256-particle count as the available FAH project 18289 input.

Only the first `getMolecules()` call is timed. Setup, ordered-output validation
and hashing are excluded. Callback counts also verify lazy computation and
cache reuse. All component counts and ordered output fingerprints match across
both DLLs and every repeat.

| Particles | Connectivity | Before seconds | After seconds | Before lifetime peak MiB | After lifetime peak MiB |
|---|---|---|---|---|---|
| 1,000,000 | water | 0.338463 | 0.042621 | 180.1 | 154.5 |
| 1,000,000 | chain | 0.185271 | 0.007402 | 171.1 | 154.5 |
| 1,000,000 | isolated | 0.326522 | 0.117198 | 177.5 | 158.4 |
| 7,950,256 | water | 2.538538 | 0.359061 | 1394.1 | 1191.1 |
| 7,950,256 | chain | 1.445829 | 0.066328 | 1316.1 | 1191.1 |
| 7,950,256 | isolated | 2.580604 | 0.883505 | 1373.3 | 1221.6 |

Windows PeakWorkingSetSize is a process-lifetime high-water mark and cannot be
reset. Raw runs retain the pre-operation working set and pre-operation lifetime
peak. In the candidate water and chain cases, the lifetime peak is set by Context
setup and does not rise during molecule discovery. Thus setup masks their actual
discovery peak. The 7,950,256-particle water-like test shows an observed whole-
process peak reduction from 1,394.1 to 1,191.1 MiB (about 203 MiB). It does not
show that the discovery operation itself consumes 1,191 MiB or that steady RAM
falls by 203 MiB. Retained RAM is essentially unchanged.

This improves the first molecule-list construction, not every cached lookup or
integration step. It is relevant to callers such as periodic wrapping and barostat setup. Static
inspection of the supplied FahCore_28.exe confirms direct calls during startup
under the messages `[ Finding non-water atoms ]` (call at image VA
`0x14001df5c`) and `[ Checking molecular weights ]` (`0x14001c950`).
These call paths share the Context cache, so this change targets the first
construction, not a repeated discovery at every checkpoint. Actual CUDA/FAH
startup time, peak RAM, checkpoint
pauses and PPD gains remain unmeasured. The platform-side graph in
`ComputeContext::findMoleculeGroups()` is unchanged; optimizing that separately
requires rebuilding and validating platform plugins. Do not add memory savings
from these measurements to earlier State-collection savings to predict one
whole-process FAH peak: the allocations occur in different phases.

## Correctness and independent review

49 executables passed: all 43 previously run serialization/State/checkpoint/
integrator suites, the existing `TestFindMolecules`, the new molecule-discovery
regression and four Reference Monte Carlo barostat suites (isotropic,
anisotropic, flexible and membrane).

The new regression compares complete ordered output with the unchanged static
algorithm on 50 deterministic random graphs, disconnected/duplicate/self edges,
a long reversed chain, constraints and nested virtual sites of all four built-in
types. It instruments custom force callbacks, verifies ordering and empty
callbacks, cached reference reuse, reinitialization and recovery after a callback
throws without publishing a partial cache.

A separate reviewer approved the design before implementation and reviewed the
final algorithm, ABI/snapshot implications, test coverage and benchmark claims.
The reviewer independently reran the new regression against both DLLs, the
candidate existing molecule test and Reference checkpoint suite. They verified
all 49 recorded outcomes, all 36 raw runs and medians, matching ordered
fingerprints, byte-identical benchmark executables and the candidate DLL hash.
The review found no correctness or ABI blocker; the setup-peak caveat above was
explicitly required for accurate reporting.

## Reproduce

Build the baseline and candidate core DLLs with the same x64 Release compiler
options. Compile the Windows benchmark once and copy the identical executable
beside both DLLs:

```powershell
cl /nologo /O2 /Ob2 /MD /EHsc /std:c++14 /DWIN32 /D_WINDOWS /I"$sourceDir/openmmapi/include" /I"$sourceDir/olla/include" /I"$sourceDir/serialization/include" /Fe:BenchmarkMoleculeDiscovery.exe "$sourceDir/devtools/benchmarks/fastxml/BenchmarkMoleculeDiscovery.cpp" "$baselineBuild/OpenMM.lib" Psapi.lib
python "$sourceDir/devtools/benchmarks/fastxml/run_molecule_discovery.py" --baseline "$baselineBuild/BenchmarkMoleculeDiscovery.exe" --candidate "$candidateBuild/BenchmarkMoleculeDiscovery.exe" --output "$resultsDir"
```

The harness requires a 6 GiB per-process committed-memory Job Object limit and
aborts if it cannot establish it. The runner applies a 180-second timeout per
process. Shared-desktop timing variation remains possible. Raw runs, medians,
compiler/build provenance and test outcomes are retained alongside this file.
