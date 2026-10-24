module PaperKernels

using LinearAlgebra, Random, Statistics, Distributions, SpecialFunctions, DelimitedFiles, Dates

# Load preserved definitions, without the original batch launch or unused imports.
function definitions_only(ex)
    ex isa Expr || return ex
    ex.head in (:using, :import) && return nothing
    ex.head == :call && ex.args[1] == :runy && return nothing
    return ex
end

for (name, source) in ((:SPMM, "spmm_entropy.jl"), (:CGPM, "cgpm_entropy.jl"),
                       (:DGPM, "dgpm_lyapunov.jl"))
    m = Module(name)
    Core.eval(m, :(using LinearAlgebra, Random, Statistics, Distributions,
                         SpecialFunctions, DelimitedFiles, Dates))
    Base.include(definitions_only, m, joinpath(@__DIR__, "..", "research", source))
    @eval const $name = $m
end

function initialize!(m, L, p)
    L isa Integer && iseven(L) && 2 <= L <= 10 || throw(ArgumentError("L must be even, 2 through 10"))
    p isa Real && isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must be in [0,1]"))
    Core.eval(m, :(nspin = $L))
    Core.eval(m, :(dimen = $(2^L)))
    Core.eval(m, :(prob = $Float64($p)))
    Core.eval(m, :(kettemp = zeros(Int, $L)))
    basis = [digits(i; base=2, pad=L) |> reverse for i in 0:2^L-1]
    Core.eval(m, :(ket = $basis))
    return nothing
end

function entropy_bits(psi, L)
    probabilities = abs2.(svdvals(reshape(psi, 2^div(L, 2), :)))
    return -sum(x * log2(x) for x in probabilities if x > 0)
end

function trajectory(model::Symbol; L=4, p=0.3, strength=0.45, steps=24, seed=2026)
    model in (:spmm, :cgpm) || throw(ArgumentError("model must be spmm or cgpm"))
    steps isa Integer && 1 <= steps <= 200 || throw(ArgumentError("steps must be 1 through 200"))
    seed isa Integer && 0 <= seed <= typemax(Int32) || throw(ArgumentError("invalid seed"))
    isfinite(strength) && strength >= 0 || throw(ArgumentError("invalid strength"))
    model == :spmm && strength > 1 && throw(ArgumentError("SPMM strength must be in [0,1]"))
    m = model == :spmm ? SPMM : CGPM
    initialize!(m, L, p)
    Random.seed!(seed)
    if model == :spmm
        Core.eval(m, :(str = $Float64($strength)))
    else
        Core.eval(m, :(mean = $Float64($strength)))
        Core.eval(m, :(Gup = Normal(mean, wid / sqrt(2))))
        Core.eval(m, :(Gdown = Normal(-mean, wid / sqrt(2))))
    end
    psi = zeros(ComplexF64, 2^L)
    psi[rand(1:length(psi))] = 1
    result = zeros(steps, 3)
    for t in 1:steps
        bonds = isodd(t) ? (1:2:L-1) : vcat(L, collect(2:2:L-2))
        for i in bonds
            psi = m.gate(psi, i)
        end
        if model == :spmm
            m.weakmsm!(psi)
        else
            psi = m.weakmsm(psi)
        end
        isfinite(norm(psi)) && isapprox(norm(psi), 1; atol=1e-10) || error("non-normalized state")
        result[t, :] = [t, entropy_bits(psi, L), norm(psi)]
    end
    return result
end

end
