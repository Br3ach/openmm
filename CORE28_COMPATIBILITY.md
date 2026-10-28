# Core28-compatible fast XML build

Branch: `fastxml-core28-8.3.1`.
Base: OpenMM tag `8.3.1`, commit `6e13f1361dcc3ffb581c7c97e3adfcc61257b813`.

This branch backports the v5 streaming State/System XML readers and moves the
large arrays through StateBuilder and Context::getState(). Existing State and
SerializationNode layouts and existing method signatures remain unchanged.
The writer retains complete SerializationNode trees, with the local allocation
improvements documented below. The v5 deferred
writer is deliberately omitted: the unchanged Core28 executable invokes its
inline serializer before a DLL-owned operation can establish a safe scope.

The incompatible master changes are upstream commits `2fbed592` (State's
hasBoxVectors member) and `6717a85c` (SerializationNode storage and property
return signatures). They are not backported. No GPU or scientific-engine
changes from newer OpenMM versions are included.

## Build

Use MSVC x64 and an out-of-tree Release build. The validated build used
MSVC 19.51, Visual Studio 18 2026 and the dynamic MSVC runtime (/MD).
For example, from the parent directory:

```powershell
cmake -S openmm-fastxml-core28-8.3.1 -B openmm-fastxml-core28-build `
  -G "Visual Studio 18 2026" -A x64 `
  -DOPENMM_BUILD_CUDA_LIB=OFF -DOPENMM_BUILD_OPENCL_LIB=OFF `
  -DOPENMM_BUILD_CPU_LIB=OFF -DOPENMM_BUILD_REFERENCE_TESTS=OFF `
  -DOPENMM_BUILD_PYTHON_WRAPPERS=OFF `
  -DOPENMM_BUILD_C_AND_FORTRAN_WRAPPERS=OFF -DBUILD_TESTING=ON
cmake --build openmm-fastxml-core28-build --config Release -j 8
ctest --test-dir openmm-fastxml-core28-build/serialization/tests `
  -C Release --output-on-failure
```

The isolated Core28 runs replaced only OpenMM.dll and retained the official
Core28 8.3.1 executable and its matched official plugin DLLs and dependencies.
No files in the installed FAH client were changed. Do not substitute plugins
built against current OpenMM master.

## Validation on 28 October 2026

- Native Release build succeeded. All 57 registered tests passed after building
  the nine optional Amoeba/Drude serialization test executables explicitly.
  This includes 36 tests in serialization/tests, the new reader regressions,
  historical force versions, unsupported-schema fallback, malformed input,
  rewind failure, complete direct proxy nodes and moved-array storage.
- A separately compiled caller using original 8.3.1 State/SerializationNode
  headers passed Context::getState(31), State return/copy/move/destruction,
  legacy property-reference lifetimes and XML round trip. Sizes remain
  State=288 bytes and SerializationNode=72 bytes on Windows x64.
- All 57 OpenMM symbols imported by the official Core28 executable are present.
- Official and candidate DLL runs completed the same generated two-particle
  WU with plain, gzip, bzip2 and encrypted gzip XML. Input integrity, startup
  scientific checks and periodic/final checks passed; all returned 100.
  XTC and force-bearing TRR bytes match. Plain/gzip/bzip2 decoded final State
  and integrator XML also match byte for byte. In this encrypted fixture,
  both official and candidate omit final State/integrator result members;
  the missing members are not claimed as an encrypted final-State comparison.
  Result packets validate and member lists match their controls.
- Temporary DLL tracing confirmed seekable Core28 input streams and successful
  fast System and State parsing, without rejection, for all four variants.
  stderr is redirected into science.log. Diagnostic source was restored and
  forcibly recompiled; final checks verify no diagnostic strings remain.
- A generated 2,000,000-step WU paused after trajectory output beyond the
  100,000-step checkpoint, exited 102, resumed and completed with 100.
  Final State, integrator, XTC and TRR bytes match the retained official
  pause/resume control from the same checksum-identical packet.

## Reader measurements

Ryzen 7 9800X3D, 100,000-particle synthetic XML; six runs per mode/role in
interleaved order, separate processes, same Release DLL and old-header caller.
OPENMM_FAST_XML=0 selects the legacy reader; =1 selects the fast reader.
Inputs are 10,500,375 bytes (State) and 2,500,276 bytes (System).
Windows process peak committed memory and working set were recorded.

| Role | Legacy median | Fast median | Legacy peak commit | Fast peak commit |
| --- | ---: | ---: | ---: | ---: |
| State | 1.038 s | 0.055 s | 219.5 MiB | 11.9 MiB |
| System | 0.205 s | 0.012 s | 65.2 MiB | 4.7 MiB |

Median peak working sets were 215.1/15.9 MiB (State) and 65.1/9.3 MiB
(System), legacy/fast respectively. These are warm filesystem-cache reader
measurements, not full Core28 startup or natural-WU performance figures.
Reported page faults include soft faults; no claim about disk paging is made.
Cold-cache behavior, decompression costs and whole-startup savings remain
unmeasured. State extraction copy reductions are functional/storage-tested;
their Core28 peak-memory benefit has not been measured separately.

## Remaining limits

Generated fixtures establish applicability, not support for every project,
plugin, checkpoint or GPU. Natural-WU startup/performance and official-to-
patched binary-checkpoint migration still need dedicated isolated tests.
The writer still retains complete particle-node trees and property maps.
Existing 8.3.1 behavior and
validation limitations outside the backported readers are unchanged.

Initial backport validation DLL SHA-256:
`d263ad380f5316689148703d87aa8ac28a46a7eb02c15becda8108d7bc99e682`.
The old-header probe, all 57 tests and all four Core28 input variants were
rerun successfully against the clean rebuild. No diagnostic strings remain.

## Writer allocation improvement (28 October 2026)

StateProxy now reads particle arrays by const reference, reserves each particle
section's known child count and reserves space for the seven possible appended
State sections in the parent. Existing children are included in the parent
reservation. This is a net four-line production increase; layouts, signatures,
section ordering and XML representation are unchanged.

All 57 tests and the original-header ABI probe passed. A new regression checks
a prepopulated parent with all seven appended sections, including integrator
parameters. Official Core28 completed the generated WU with exit 100; final
State, integrator, XTC and TRR bytes match the official control.

A separate executable compiled against original 8.3.1 headers serialized
100,000-particle States through the baseline and optimized DLLs. Six runs per
case/build were interleaved in alternating order, using separate processes.
Time includes XML formatting and file output; memory includes State setup.
XML files matched byte for byte for each case. Median measurements:

| State contents | Baseline time | Optimized time | Baseline peak commit | Optimized peak commit |
| --- | ---: | ---: | ---: | ---: |
| Positions | 0.322 s | 0.194 s | 112.6 MiB | 61.2 MiB |
| Positions, velocities, forces, parameters, energy | 0.983 s | 0.574 s | 231.3 MiB | 176.6 MiB |
| All above plus integrator parameters | 1.104 s | 0.572 s | 343.6 MiB | 176.6 MiB |

These synthetic serialization measurements are not predictions for a natural
7.97-million-particle WU. They measure the combined reference/reservation
change, not the separate contribution of each reservation.

Writer-optimized DLL SHA-256:
`f626478061e487fc3abbc7829bb0bdfb1441d332f80896a7657d4dd21060e91e`.
