# Measured results

Baseline: original fastxml-clean-v4 commit `655c38e9179d984f095f32d811d3eceed7c69534`. Candidate: the implementation and source hashes recorded in `validation.json`.

Medians of three fresh-process runs. RAM is peak resident memory, in MiB. Both harnesses were compiled separately against their matching headers.

| Fixture | Operation | Baseline peak MiB | Candidate peak MiB | Baseline seconds | Candidate seconds |
|---|---|---:|---:|---:|---:|
| system | read | 73.3 | 73.3 | 0.541 | 0.568 |
| system | write | 911.7 | 76.6 | 2.850 | 0.555 |
| system-fallback | read | 931.7 | 73.9 | 4.806 | 0.594 |
| system-fallback | write | 911.6 | 76.5 | 2.639 | 0.535 |
| state-pv | read | 50.4 | 50.5 | 0.574 | 0.604 |
| state-pv | write | 53.5 | 53.4 | 1.321 | 1.399 |

The System fixtures contain 500,000 particles; the State contains 1,000,000 particles with positions and velocities. The late-fallback fixture appends a small CustomBondForce. See the parent README for methodology and limitations.

39/39 serialization tests passed. All three checked large XML outputs were byte-identical. Unsupported force subtrees may still be large; scanner-incompatible XML retains full fallback. Caller recompilation is required for the new System writer. These results do not measure simulation throughput or memory compression.

Native Windows results and real FAH XML validation: [windows/README.md](windows/README.md).
