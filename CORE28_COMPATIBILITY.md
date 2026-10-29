# Core28-compatible fast XML build

Branch: `fastxml-core28-8.3.1`.
Base: OpenMM tag `8.3.1`, commit `6e13f1361dcc3ffb581c7c97e3adfcc61257b813`.

This branch backports the v5 streaming State/System XML readers and moves the
large arrays through StateBuilder and Context::getState(). Existing State and
SerializationNode layouts and existing method signatures remain unchanged.
The ordinary writer retains complete SerializationNode trees. A narrowly
guarded adapter streams State arrays for the verified Windows Core28 8.3.1
executable only; its compatibility boundary and measurements are documented
below. No public layout or method signature is changed.

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
plugin, checkpoint or GPU. Whole-Core28 startup/performance and official-to-
patched binary-checkpoint migration still need dedicated isolated tests.
Outside the verified Core28 State adapter, the writer retains complete
particle-node trees and property maps.
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

## Legacy State deserialization improvement (28 October 2026)

The legacy StateProxy reader reserves each particle array from the child count,
records its size before moving it into StateBuilder, moves the parameter map,
and returns the constructed State through takeState(). Existing count checks,
duplicate-section behavior and the const input tree are preserved. Production
code grows by two lines. Successful fast-reader calls are unaffected.

All 57 tests and the original-header ABI probe passed. Regressions cover unequal
array lengths, unequal duplicate sections, empty-array flags, input-node reuse,
IntegratorParameters fallback and cloning. Official Core28 completed the
isolated generated WU with OPENMM_FAST_XML=0 and exit 100; final scientific
and trajectory members are byte-identical to the official control.

Six interleaved baseline/candidate runs per case used the same 100,000-particle
inputs and MSVC build settings. For direct proxy deserialization from an
identical prebuilt node tree, the full State median was 0.0697 -> 0.0644 s;
process peak committed memory was 192.7 -> 183.6 MiB. With integrator parameters,
medians were 0.0690 -> 0.0666 s with the same approximate memory reduction.
These process peaks include construction of the retained source tree.

Full legacy XML deserialization (OPENMM_FAST_XML=0) showed no material gain:
the full State medians were 1.057 -> 1.086 s, with approximately 228 MiB peak
commit in both builds. With integrator parameters both medians were 1.233 s
and approximately 343 MiB. XML tree construction dominates these cases; no
whole-startup speedup is claimed. The local change is retained for its measured
direct-proxy memory reduction and removal of unnecessary copies without a new
abstraction or altered validation behavior.

Legacy-reader-optimized DLL SHA-256:
`086c73eac15f447caf8752573b06aa97ea8d48676ad1f8ada62901365306f7bf`.

## Buffered tag-body appends (28 October 2026)

FastXmlInput now scans contiguous buffer spans and appends each span to the
existing tag string. It replaces the old per-character tag-body loop; prefix,
whitespace, BOM, comment and declaration handling remain unchanged. Quote
state survives refills. Scanning stops at the first unquoted '>' or the first
byte beyond the existing tag-size limit, preserving consumption and rewind.
The change adds one internal helper and is a net 17 production lines.

All 57 tests and the original-header ABI probe pass. The new boundary suite
also passes against the unchanged baseline DLL. It exercises both quote forms,
quoted '>', opening/closing quotes and delimiters at 1 MiB boundaries +/-1,
short refills of 1/2/3/7/31 bytes, BOM/comments/declarations, exact tag-size
limits, truncation, and trailing newline/data. Its failed-rewind guard ensures
successful cases use the fast reader rather than silently falling back.
Existing tests compare canonical serialized State/System contents with the
stock parser, including numeric/force properties and unsupported-schema fallback.

Official Core28 completed the isolated generated WU for plain, gzip, bzip2
and encrypted inputs, all with exit 100. Scientific/trajectory result members
match the controls byte for byte where present; encrypted controls omit the
final XML members as documented above. No installed FAH files were changed.

Six same-session, interleaved baseline/candidate pairs parsed the retained
1,927,740,675-byte natural System XML (SHA-256 abda690dfe141b24284d7a9288840804eedf9603bcce5a5af607c3c1db481ca0).
The helper, text-mode ifstream, input, MSVC Release /O2 /Ob2 settings and machine
were identical; no context or simulation was created. The direct fast-reader
entry avoids conflating these measurements with a legacy fallback. All runs
returned 7,950,180 particles, 5,955,581 constraints, six forces and 11,940,050
nonbonded exceptions. Counts alone do not establish complete natural-System
numerical equality; semantic/output checks above use the isolated test fixtures.

| Natural System parsing | Baseline median | Optimized median |
| --- | ---: | ---: |
| Wall time | 17.377 s | 14.183 s |
| CPU time | 17.070 s | 14.047 s |
| Peak committed memory | 1.5933 GiB | 1.5933 GiB |

Wall time decreased 18.4%, CPU time 17.7%; every pair favored the candidate.
These same-session controls are the comparison, not earlier profiling timings
from another session. No cache flush or CPU-frequency lock was applied; this
is not a cold-cache or complete Core28 startup speedup claim.

A separate six-pair 100,000-particle State check gave median wall times
0.06437 -> 0.05409 s with effectively unchanged peak memory. Independent
review found no actionable findings and supported retaining the localized
change after these validation results.

Tag-span DLL SHA-256:
`7b134f840aab850dac640bd8fd761960391d20d07812be17f290730b929fa48d`.

## In-place child-node construction (28 October 2026)

SerializationNode::createChildNode() now uses emplace_back() instead of
push_back(SerializationNode()), preserving subsequent name assignment and the
returned reference. This is one replacement line with zero production growth.
All 57 tests and the original-header ABI probe passed. Existing writer tests
produced byte-identical XML for positions-only, full State and full State with
integrator parameters in six interleaved baseline/candidate pairs each.
Median writer times were respectively 0.23684 -> 0.23629 s, 0.67479 -> 0.69877 s
and 0.70332 -> 0.68533 s. Memory was effectively unchanged. These mixed results
do not establish a serialization speedup; the change is retained as a local
direct-construction simplification, not a measured performance improvement.
Independent review found no actionable findings.

DLL SHA-256:
`c78d906a236db1a7c17131ee438e86bf69808ce871f35a7c7c6e7be3aa9d3b9d`.

## XML escaping experiments (28 October 2026)

The retained change replaces the five-character escape map with a switch,
removing the mutable static map and its manual initialization flag (net three
fewer production lines). Hexadecimal-reference pass-through, control-byte
formatting and raw high-byte behavior are preserved, including existing
handling of malformed references. This is a small stateless simplification;
the measurements below do not establish a writer speedup.

A separately evaluated string-reuse candidate cleared name/value per attribute
inside a block ending before child recursion. It was rejected because its
added scope and lifetime bookkeeping produced no consistent measured benefit.
The retained code still uses the original per-attribute temporary strings.
No direct-stream writer or packed-node representation was introduced.

Each phase used six interleaved baseline/candidate pairs per 100,000-particle
writer case. Time includes XML formatting and file output. Median seconds:

| Case | Original map | Switch (phase 1) | Switch (phase 2 control) | Switch + reuse |
| --- | ---: | ---: | ---: | ---: |
| Positions | 0.23124 | 0.24232 | 0.27017 | 0.26256 |
| Full State | 0.72147 | 0.73353 | 0.68981 | 0.67946 |
| Full State + integrator parameters | 0.74928 | 0.74586 | 0.67482 | 0.67990 |

Only compare columns within their phase; the separate phases had independent
controls. Differences were mixed/small, and peak memory was effectively
unchanged. The large-State XML files were byte-identical in both phases.

The byte oracle serialized 781 inputs as both values and property names
(1,562 documents), covering all bytes 0..255, empty/long strings, UTF-8,
named/decimal/uppercase references, and malformed lowercase hexadecimal
prefixes. Original, switch and reuse builds produced identical oracle bytes
(SHA-256 c72ec793d9685663e456ba21d1de3e66ded0743889544980bcc8e89e313a0668).
Regression tests also explicitly check the five escapes, all control bytes
including embedded NUL, high bytes and the legacy reference quirks. These
are byte comparisons; invalid XML output is not required to round-trip.

All 57 tests and the original-header ABI probe passed on the retained build.
The new test suite also passed against the unchanged baseline DLL. Official
Core28 completed the generated plain-input WU with matching final scientific
and trajectory member bytes. Independent review supported retaining only
the stateless switch and dropping reuse, without a performance claim.

Retained switch-only DLL SHA-256:
`2c6ee7849720f93267c1632be88e728b4c38b2662aa6648685b6f64cf77f6ff5`.

## Three-particle virtual-site parsing (28 October 2026)

The two retained natural Systems contain 1,984,130 and 1,955,136
ThreeParticleAverageSite elements respectively, and no other virtual-site
types. The first input is 1,927,740,675 bytes, SHA-256
`abda690dfe141b24284d7a9288840804eedf9603bcce5a5af607c3c1db481ca0`.

Temporary timing of the complete virtual-site function attributed 1.898–2.007
seconds of a 13.763–13.859-second parse to this path. The empty timer-pair
calibration loop took 0.078 seconds for the same call count; summed empty
clock gaps were 0.039 seconds. These are inclusive diagnostic estimates,
including construction, allocation and leaf consumption. Timer overhead and
changed inlining prevent precise attribution to attribute lookups alone.
The diagnostic build's whole parse, excluding calibration, was comparable
to the unchanged build; no diagnostic code is retained.

The retained candidate replaces six complete attribute scans with one pass
using the existing FastXmlAttributes iterator, duplicate checks and numeric
converters. Only ThreeParticleAverageSite changes, at 14 net production lines.
Required fields, ignored unknown attributes, construction and leaf handling
are preserved. No class layouts, public interfaces or other site types change.

Six sequential warm-cache baseline/candidate pairs used the same MSVC Release
settings and natural System. Order alternated between pairs. Seconds:

| Pair | First | Baseline wall | Candidate wall |
| --- | --- | ---: | ---: |
| 0 | Baseline | 14.033311 | 11.871398 |
| 1 | Candidate | 13.798692 | 12.078169 |
| 2 | Baseline | 14.003782 | 11.990353 |
| 3 | Candidate | 14.642723 | 12.167744 |
| 4 | Baseline | 14.052084 | 12.331912 |
| 5 | Candidate | 14.050586 | 12.339392 |

Median wall time fell from 14.041949 to 12.122957 seconds (13.7%); median
CPU time fell from 13.984375 to 12.046875 seconds (13.9%). All six pairs
improved. Median peak private commit was effectively unchanged, 1,710,802,944
versus 1,710,718,976 bytes; median peak working set differed by only 38 KiB.
These measure System parsing, not complete Core28 startup or GPU throughput.
No cold-cache claim is made; caches were not flushed and live folding was
not interrupted.

Outside the timed interval, every virtual site's particle ordinal, three
parent indices and three double-precision weights were hashed in host byte
order with FNV-1a64. All twelve parses produced `cce611c4f214eca5`, with
1,984,130 sites, 7,950,180 particles, 5,955,581 constraints and six forces.
This digest is a regression comparison, not a cryptographic integrity check.

One additional pair on the second retained System parsed in 14.216050 versus
12.203553 seconds. Both builds returned 7,863,224 particles, 5,868,599
constraints, six forces and 1,955,136 sites, with site digest
`7f9db4271e813845`. This is a compatibility check and a single timing pair,
separate from the six-pair benchmark above.

All 57 existing tests, the Reference virtual-site scientific test and the
original-header ABI runtime probe passed. Focused tests compare stock and
fast-parsed System XML and prohibit silent fallback for valid cases. They
cover attribute order, quotes, unknown and duplicate unknown attributes,
explicit closing tags, seven-byte refills, missing/duplicate required fields,
integer limits, invalid numbers and malformed attributes. The same focused
test executable also passed against the original DLL.

Validated candidate DLL SHA-256:
`e7c72bd584eb7e97458eb26c519088e63ae2841a9523c84452b6f9b515bc0e89`.

## Proxy and text-stream compatibility fixes (28 October 2026)

Fast State/System dispatch now checks the exact dynamic type of the registered
proxy against the existing built-in proxy class. The System reader performs
the same check for each supported force it encounters. Overrides, including
derived proxy classes, use the existing checked fallback and registered
deserializer. An override of an unused force does not disable fast parsing.
The registry itself and all public interfaces/layouts are unchanged.

Read-ahead restoration previously subtracted an unread character count from a
file position. Windows text-mode CRLF conversion made that incorrect: a tail
of `\r\nTAIL\r\nSECOND\r\n` could become `AIL\nSECOND\n`. The reader now
saves the opaque stream position before each refill and restores it before
replaying the consumed decoded characters through the existing buffer.
Replay handles positive short reads and rereads at most one 1 MiB buffer;
failed saves, seeks or reads retain the existing error/fallback handling.
There is no extra XML-sized buffer and no text/binary-mode heuristic.

The two fixes add 21 net production lines. Tests now observe the legacy
reader's end-position query instead of assuming that every fourth positional
seek is fallback: successful fast restoration itself legitimately seeks.
Coverage includes State, System and all six supported force proxy overrides;
an unused override; Windows text/binary files with LF/CRLF XML, trailing data
and buffer-spanning documents; short replay reads, replay EOF/read exceptions
and a failed finish seek. All 57 tests, both Reference nonbonded/virtual-site
scientific tests and the original-header ABI runtime probe passed.

Four sequential warm-cache pairs compared the preceding commit `c95e8f2d`
against the fixed build, alternating order on the same natural System:

| Pair | First | Before wall (s) | After wall (s) |
| --- | --- | ---: | ---: |
| 0 | Before | 11.894577 | 11.828031 |
| 1 | After | 12.010728 | 11.541027 |
| 2 | Before | 12.078664 | 11.997217 |
| 3 | After | 11.837520 | 12.036507 |

Median wall time was 11.952653 versus 11.912624 seconds; CPU time was
11.914063 versus 11.812500 seconds. Peak private commit was effectively
unchanged (1,710,761,984 versus 1,710,798,848 bytes). This establishes no
meaningful regression in this System-parsing workload, not a speedup claim.
All eight runs retained the full virtual-site digest `cce611c4f214eca5`.

Official Core28 completed fresh isolated generated WUs using plain, gzip,
bzip2 and encrypted gzip XML, all with exit 100 and valid result packets.
Returned State/Integrator XML and XTC/TRR members matched the retained official
outputs byte for byte after decompression. The encrypted case retained the
official omission of State/Integrator members and matched its trajectory
members. Logs contain run-specific values and are not asserted byte-identical.
These are tiny generated-WU checks, not natural-WU completion or throughput
measurements. Installed FAH files and live WUs were not changed.

Validated DLL SHA-256:
`5b57067309e34d9c5ab98cca0e489e35159b7308df721b0146ac3f7e32e15b58`.


## Fixed-binary Core28 State adapter (29 October 2026)

This adapter is restricted to MSVC x64 and the official Core28 8.3.1 executable
with SHA-256
`9e27f4493c83a31e8fb0403248629909a67f8fc9c3a7ac6a41afa6dc93133436`.
It verifies the on-disk image, the loaded serializer wrapper and both call
sites. A verification failure, including an unavailable or overlong executable
path, retains the ordinary complete-tree path before any deferral. The existing
`OPENMM_FAST_XML=0` switch disables the adapter too.

Only the built-in State proxy, its fresh canonical root, and the observed
nonempty matching position/velocity/force arrays with State flags `0x1f` are
eligible. The verified wrapper writes synchronously while that State remains
alive. A thread-local handoff is consumed before output; allocation failure in
the intervening type-property insertion cancels it without allocating.
Replacement of the registered State proxy permanently disables future deferrals
for that process, even if the original is restored. An existing handoff remains
valid. This protects custom proxies that tail-call the built-in proxy. Ordinary
direct proxy use still produces complete nodes; a similarly named user property
alone never activates streaming.

The adapter reuses XML escaping and double formatting, writes at most 4 MiB of
buffered output through the original ostream, and preserves stream-state and
exception-mask behavior. Core28 retains compression, checkpoint ordering,
scientific checks and file ownership. SerializationNode constructors,
destructors, copying and layout are unchanged. Windows SHA verification uses
the operating system's `bcrypt.dll`; no additional library needs shipping.
The adapter does not promise recovery from asynchronous faults or memory
corruption. It is tied to this exact executable, not a general lazy-node API.

Validation of optimized DLL SHA-256
`358416fa766200c57a6d3083efa9d7f9ac49769fca5f2c1ef0b4553112417f15`:

- All 57 non-Reference native tests and the Reference nonbonded/virtual-site
  tests pass. Original-header ABI and API oracle probes pass, with unchanged
  exports and all imports of the official executable and four plugins present.
- Official Core28 with the final DLL completes a generated WU and a clean
  pause (102), resume and completion (100). Decoded final State, integrator,
  XTC and TRR match the retained controls byte for byte. DLL-only debugger
  tracing confirms that the ordinary matched host actually registers a handoff.
- Active custom proxy and restore-before-tail-call cases complete and match
  their baseline outputs. The former reproduced a failure before the sticky
  opt-out was added. A changed-executable negative control retains the tree.
- Plain, gzip, bzip2 and encrypted generated inputs passed again on the final
  optimized DLL and match retained official controls. Encrypted fixtures omit
  final State/integrator members in both controls, so their State output parity
  is not established by that test.
- A separately linked diagnostic build of the same source exercises handoff
  cancellation/reuse/nesting, thread isolation, registration/setter allocation
  failures, reentrant new_handler behavior, byte formatting and ostream failures.
  Diagnostic builds are not used for performance measurements.

Three alternating baseline/candidate pairs on a generated 250,000-particle WU
measured median total run time of 4.019 -> 2.580 seconds, peak private memory
1483.21 -> 1033.30 MiB and peak working set 796.10 -> 363.53 MiB. These pairs
precede only the proxy-history safeguard. A final-build confirmation pair gave
3.931 -> 2.692 seconds, 1482.03 -> 1031.77 MiB private and 796.11 -> 362.18 MiB
working set. Both final result XML members match byte for byte.

These are warm-cache, isolated full-run measurements, including startup,
validation, simulation and packaging, not timings of XML alone. The baseline
is this branch's preceding DLL (`0e0f7ea6...`), not stock OpenMM. No cache eviction
or live-installation modifications were used. Natural-WU peak savings, cold
cache behavior and broad project compatibility remain unmeasured; the generated
results must not be extrapolated into an 8-million-particle memory guarantee.
