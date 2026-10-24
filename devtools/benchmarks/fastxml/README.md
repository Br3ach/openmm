# XML peak-memory changes and benchmark

The `fastxml-clean-v4` baseline (`655c38e9179d984f095f32d811d3eceed7c69534`)
streams supported System input and large State input/output, but still builds
the entire System serialization tree on output. One unsupported force on input
also restarts the entire System through the generic reader.

This change targets those two remaining peaks. It does not introduce memory
compression, change simulation kernels, or optimize binary checkpoints.

## Design

### Stream top-level System output

The `XmlSerializer::serialize()` template calls an internal helper before
building a proxy tree. The helper accepts only an exact `System` with its exact
built-in `SystemProxy`, canonical root name `System`, at least 100,000 particles
or constraints, and the classic global C++ locale. All other calls retain the
original proxy path. Both `OPENMM_FAST_XML=0` and
`OPENMM_FAST_SYSTEM_WRITE=0` disable this writer.

The writer uses a checked 4 MiB output buffer. It emits particles, virtual sites,
constraints, and the large arrays of the exact built-in NonbondedForce,
HarmonicBondForce, HarmonicAngleForce, and PeriodicTorsionForce proxies directly.
It retains stock property ordering, floating-point formatting, escaping and
empty-element formatting. Other forces, including plugins and replaced proxies,
are serialized one at a time through the registered proxy. No output retry is
attempted after a write failure.

An independent design review rejected a borrowed-pointer handoff from
`SystemProxy::serialize()`: proxy nodes may be nested, moved, retained, or read
after the original System has been destroyed. The template hook preserves those
snapshot semantics. `SystemProxy` and its generic tree contract are unchanged.

**Callers must be recompiled against the updated `XmlSerializer.h` to benefit
from System writing.** Existing binaries remain usable with the new library,
but their System-writing template instantiation retains the original path.
The added helper changes no existing class layout or virtual table. Newly built
callers require the updated library providing the helper symbol.

### Limit input fallback to one force

The reader tracks logical stream positions independently of its 1 MiB read-ahead
buffer. Before parsing a force, it records the force's opening position. If the
direct parser rejects that force, or its registered proxy has been replaced,
it rewinds only to that force. A matching-element stack collects its XML into an
isolated stream and the existing generic reader invokes the registered proxy.
The surrounding System continues streaming afterward.

The isolated stream is necessary because the legacy generic reader seeks to
absolute stream offset zero. Existing irrXML handles attribute decoding, and the
collector discards comments and processing instructions, as generic nodes do.
A completed generic proxy failure propagates its original exception; it does
not invoke the proxy a second time by retrying the entire System. Exact root
proxy checks also preserve custom State/System deserialization.

Whole-document fallback remains for XML constructs the fast scanner cannot
represent (including text/CDATA), or when a safe local rewind cannot be made.
State IntegratorParameters are unchanged. Generic fallback memory is bounded
by the unsupported force subtree plus the live System, not by a constant amount
independent of force size. Caller-retained XML buffers also remain outside this
optimization.

## Validation

The serialization tests cover canonical byte equality; all direct force writers;
virtual sites, offsets, escaped names and empty sections; custom root/force
proxies; clone and standalone proxy snapshots; nested large Systems; output
failures; nonclassic locales; mixed force order; known-force schema rejection;
fragment and rewind boundaries across 1 MiB buffers; nonzero document offsets;
preserved stream suffixes; denied local rewinds; and single invocation of a
throwing custom proxy.

See `results/` for the measured baseline/candidate runs, source-file hashes and
validation outcomes. These are real executions on synthetic inputs. They do not
reproduce a specific Folding@home work unit, CUDA Context setup, Windows memory
management or OS memory compression.

## Reproduce the memory benchmark (Linux)

Build the original branch commit and the candidate with the same Release
configuration. For a CPU-only serialization benchmark, disable CUDA, OpenCL,
HIP, CPU, Amoeba, RPMD, Drude and PME plugins, C/Fortran/Python wrappers and
examples. Build the `OpenMM` target. Compile this harness **separately against
each checkout's headers** and link it to that checkout's library:

```sh
g++ -O3 -DNDEBUG -std=c++17 BenchmarkXmlMemory.cpp \
  -I"$source_dir/openmmapi/include" -I"$source_dir/olla/include" \
  -I"$source_dir/serialization/include" -L"$build_dir" \
  -Wl,-rpath,"$build_dir" -lOpenMM -o "$benchmark_executable"
```

Use absolute paths for the source, build and output executable variables. Then:

```sh
python3 run_benchmark.py --baseline /absolute/path/bench-baseline \
  --candidate /absolute/path/bench-candidate --output /absolute/path/results
```

Requirements: Linux `/proc/self/status` and `/proc/self/clear_refs`, glibc
`malloc_trim`, GCC-compatible C++17 compiler, Python 3.11+, and approximately
3 GiB available RAM. The runner enforces a 6 GiB per-process virtual-address
space limit and a 180-second per-operation timeout. Each version/operation runs
three times in fresh processes, alternating variant order. Files are warmed
before reading; writes go to `/dev/null`. These timings exclude cold-disk and
durable-storage costs.

The System fixture has 500,000 particles in water-like groups of three, a
NonbondedForce, and three constraints/excluded pairs per complete group. The
fallback fixture appends a small CustomBondForce. The State control contains
positions and velocities for one million particles. These are serializer
workloads, not validated molecular simulations.

Write inputs are constructed before measurement, released setup allocations
are trimmed, and Linux VmHWM is reset before timing. Operation peak RSS includes
live inputs and library overhead. Baseline and retained RSS, plus lifetime
`getrusage` peak, are recorded separately. Retained RSS can include allocator
caches; vector capacity is not equivalent to resident RAM.

Every run validates counts and hashes the synthetic numerical data. Candidate
XML output is also checked byte-for-byte against the baseline fixtures after
timing. Source and outputs are not included in elapsed time. Shared-host timing
variation remains possible; raw runs, medians and ranges are retained.
