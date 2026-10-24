# Paper Results and Source Map

[Overview](../README.md) | [Code interfaces](CODE.md) | [Reproducibility](REPRODUCIBILITY.md)

Figure numbers refer to [arXiv:2404.02968v2](https://arxiv.org/abs/2404.02968v2).
The images below are published results, not newly simulated data.

## Phase Diagrams: Figure 2

![Figure 2](figures/fig02.png)

Measurement strength and measurement rate are varied independently. The color
maps show entropy variance; the overlaid points use ancilla probes to locate
the boundary between volume-law (VL) and area-law (AL) regimes.

[Original PDF](figures/fig02.pdf).
The recovered [CGPM](../research/cgpm_entropy.jl) and
[SPMM](../research/spmm_entropy.jl) templates compute entropy-variance observables.
The [quickstart](../examples/quickstart.jl) calls their local dynamics kernels
but records a short entropy trace instead of scanning a phase diagram.

## Ancilla Probe: Figure 3

![Figure 3](figures/fig03.png)

The ancillary-qubit entropy provides a transition diagnostic. Finite-size
collapse and time dependence are used to investigate critical behavior.

[Original PDF](figures/fig03.pdf).
Recovered related code: [SPMM ancilla driver](../research/spmm_ancilla.jl).
[CGPM bipartite](../research/cgpm_bipartite.jl) and
[two-ancilla correlation](../research/cgpm_ancilla_correlations.jl) templates
are also preserved, but these are different observables, not interchangeable
drivers for this entire figure. Exact production-script-to-panel provenance
has not been established for every panel.

## Gaussian Transfer Matrix: Figure 5

![Figure 5](figures/fig05.png)

The measurement record and two-vector evolution probe free-energy scaling,
Lyapunov gaps, and their fluctuations. Discretizing the pointer outcomes is
important for the measurement-record entropy used here.

[Original PDF](figures/fig05.pdf).
Code: [DGPM Lyapunov template](../research/dgpm_lyapunov.jl).
It writes the leading accumulated cost, a two-vector gap, and the squared gap.
The [tests](../tests/runtests.jl) exercise this driver at L=4.
The archive's default p=0.194 is preserved; the figure caption specifies p=0.19.
The [strength-four variant](../research/dgpm_str4_lyapunov.jl) is a related
archived Gaussian calculation, not the strength-one production configuration.

## Softened Transfer Matrix: Figure 6

![Figure 6](figures/fig06.png)

This is the corresponding published transfer-matrix analysis for softened
projective measurements.

[Original PDF](figures/fig06.pdf).
The [SPMM free-energy template](../research/spmm_free_energy.jl) is recovered.
The complete original two-vector SPMM campaign and postprocessing are not yet
identified in this repository; do not interpret the presence of this figure
as a claim that the quickstart reproduces its multifractal panels.

## What You Can Reproduce Now

| Item | Available here |
| --- | --- |
| Small CGPM/SPMM trajectories and entropy trace | Yes: [quickstart](../examples/quickstart.jl), [saved outputs](../examples/expected/) |
| Small DGPM two-vector run | Yes: [test](../tests/runtests.jl) |
| Original archived simulation definitions | Eight source files with [hashes and origins](source_manifest.json) |
| Published figure previews and original figure assets | Figures 2, 3, 5, 6, with [attribution](figures/README.md) |
| Full ensemble data and all final fitting pipelines | Not included |

No critical exponent should be fitted from the tiny demonstration traces.
