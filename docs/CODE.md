# Simulation Interfaces

## Supported Example API

```julia
include("examples/kernels.jl")
rows = PaperKernels.trajectory(:spmm; L=4, p=0.28, strength=0.45,
                               steps=24, seed=2026)
```

`model` is `:spmm` or `:cgpm`; the result has columns `layer`,
`halfcut_entropy_bits`, and `state_norm`. One layer means one alternating
brickwork gate layer followed by measurements. Spatial boundaries are periodic.
CGPM strength means pointer separation / width; SPMM strength is Lambda.
The former is nonnegative and the latter is restricted to [0,1].

The simulator keeps a full complex statevector. The example uses an initial
computational-basis product state. Measurements sample their outcomes from the
Born probabilities and normalize the conditional state; this is not forced
postselection. Entropies use base-2 logarithms; the DGPM measurement-record costs
use natural logarithms.

The wrapper validates even L in 2:10, probability in [0,1], layer count in 1:200,
and a nonnegative 32-bit seed. It resets the default RNG and module state for
each call. Calls must be serial; separate processes can isolate independent runs.

## Archived Sources

| File | Role |
| --- | --- |
| [cgpm_entropy.jl](../research/cgpm_entropy.jl) | Continuous Gaussian-pointer evolution and entropy variance |
| [spmm_entropy.jl](../research/spmm_entropy.jl) | Binary softened-measurement evolution and entropy variance |
| [dgpm_lyapunov.jl](../research/dgpm_lyapunov.jl) | Discrete-pointer two-vector evolution and QR orthogonalization |
| [dgpm_str4_lyapunov.jl](../research/dgpm_str4_lyapunov.jl) | Gaussian strength-four variant |
| [spmm_free_energy.jl](../research/spmm_free_energy.jl) | Softened-measurement record cost |
| [spmm_ancilla.jl](../research/spmm_ancilla.jl) | Ancilla entropy scan |
| [cgpm_bipartite.jl](../research/cgpm_bipartite.jl) | Bipartite entropies |
| [cgpm_ancilla_correlations.jl](../research/cgpm_ancilla_correlations.jl) | Two-ancilla correlations |

These are research templates, not a polished package API. Their original global
variables, historical output paths, and batch calls are preserved. **Do not run
the archived files directly as the quickstart.** Some defaults start larger
jobs or write to historical locations.

`examples/kernels.jl` uses Julia's parsed-expression `include` hook to suppress
top-level `runy(...)` calls and original import declarations. It supplies the
runtime dependencies needed by the selected kernels explicitly. It does not
alter function bodies or evaluate generated strings. Unused CSV, DataFrames,
StatsBase, and Bootstrap imports are not required by these selected examples.
The example schedules short trajectories and computes entropy with an independent
Schmidt decomposition. Tests compare it with the original entropy routine.

Sources are hashed in [source_manifest.json](source_manifest.json). Their presence
does not certify every old helper or hard-coded batch configuration.
