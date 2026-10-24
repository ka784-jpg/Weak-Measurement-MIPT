using Test, LinearAlgebra, Random, DelimitedFiles
include("../examples/quickstart.jl")

@testset "Original gate and measurement kernels" begin
    for m in (PaperKernels.SPMM, PaperKernels.CGPM, PaperKernels.DGPM)
        for seed in 1:5
            Random.seed!(seed)
            U = m.rangate()
            @test U' * U ≈ I atol=1e-12
            reshuffled = reshape(permutedims(reshape(U, 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
            @test reshuffled' * reshuffled ≈ I atol=1e-12
        end
    end
    for model in (:spmm, :cgpm)
        a = PaperKernels.trajectory(model; L=4, steps=12, seed=9)
        @test a == PaperKernels.trajectory(model; L=4, steps=12, seed=9)
        @test all(isfinite, a)
        @test all(x -> -1e-12 <= x <= 2 + 1e-12, a[:, 2])
        @test all(x -> isapprox(x, 1; atol=1e-12), a[:, 3])
        m = model == :spmm ? PaperKernels.SPMM : PaperKernels.CGPM
        PaperKernels.trajectory(model; L=4, p=1, strength=0, steps=1)
        psi = normalize(randn(ComplexF64, 16))
        before = copy(psi)
        after = model == :spmm ? (m.weakmsm!(psi); psi) : m.weakmsm(psi)
        @test after ≈ before atol=1e-12
    end
    rows = PaperKernels.trajectory(:spmm; L=4, p=1, strength=1, steps=8)
    @test maximum(abs, rows[:, 2]) < 1e-10
    psi = normalize(randn(ComplexF64, 16))
    PaperKernels.initialize!(PaperKernels.SPMM, 4, 0.0)
    @test PaperKernels.SPMM.reyni(psi, 1, 2) ≈ PaperKernels.entropy_bits(psi, 4)
end

@testset "Dense gate reference" begin
    m = PaperKernels.SPMM
    PaperKernels.initialize!(m, 4, 0.3)
    Random.seed!(17)
    psi = normalize(randn(ComplexF64, 16))
    for i in (1, 2, 3)
        Random.seed!(22)
        U = m.rangate()
        dense = kron(Matrix{ComplexF64}(I, 2^(4-i-1), 2^(4-i-1)), U,
                     Matrix{ComplexF64}(I, 2^(i-1), 2^(i-1)))
        Random.seed!(22)
        @test m.gate(copy(psi), i) ≈ dense * psi atol=1e-12
    end
end

@testset "DGPM archived driver smoke test" begin
    mktempdir() do dir
        Random.seed!(51)
        m = PaperKernels.DGPM
        m.runy(4, 0.194, 1, joinpath(dir, "sample"))
        values = readdlm(joinpath(dir, "sampleiterat=1_FET"), ',', ComplexF64)
        @test size(values) == (112, 3)
        @test all(isfinite, values)
        @test maximum(abs, imag.(values)) < 1e-10
        @test m.state' * m.state ≈ I atol=1e-10
    end
end

@testset "Input validation and safe outputs" begin
    @test_throws ArgumentError PaperKernels.trajectory(:unknown)
    for kwargs in ((; L=3), (; p=NaN), (; p=-0.1), (; strength=2), (; steps=0), (; seed=-1))
        @test_throws ArgumentError PaperKernels.trajectory(:spmm; kwargs...)
    end
    mktempdir() do dir
        out = joinpath(dir, "demo")
        main([out])
        @test isfile(joinpath(out, "spmm.csv"))
        @test isfile(joinpath(out, "cgpm.csv"))
        saved = read(joinpath(out, "spmm.csv"))
        @test_throws ArgumentError main([out])
        @test saved == read(joinpath(out, "spmm.csv"))
    end
end
