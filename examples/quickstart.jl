using DelimitedFiles, TOML, SHA
include("kernels.jl")
using .PaperKernels

function main(args)
    length(args) == 1 || throw(ArgumentError("usage: julia --project=. examples/quickstart.jl OUTPUT_DIRECTORY"))
    output = abspath(only(args))
    ispath(output) && throw(ArgumentError("output already exists: $output"))
    parent = dirname(output)
    mkpath(parent)
    stage = mktempdir(parent)
    try
        for (model, strength, p) in ((:spmm, 0.45, 0.28), (:cgpm, 1.0, 0.19))
            rows = PaperKernels.trajectory(model; L=4, p, strength, steps=24, seed=2026)
            open(joinpath(stage, "$model.csv"), "w") do io
                println(io, "layer,halfcut_entropy_bits,state_norm")
                writedlm(io, rows, ',')
            end
        end
        metadata = Dict("L" => 4, "layers" => 24, "seed" => 2026,
                        "models" => Dict("cgpm" => Dict("strength" => 1.0, "p" => 0.19),
                                         "spmm" => Dict("strength" => 0.45, "p" => 0.28)),
                        "julia_version" => string(VERSION),
                        "description" => "Two single trajectories, not published ensembles",
                        "sources" => Dict(file => bytes2hex(sha256(read(joinpath(@__DIR__, "..", "research", file))))
                                          for file in ("spmm_entropy.jl", "cgpm_entropy.jl")))
        open(joinpath(stage, "metadata.toml"), "w") do io
            TOML.print(io, metadata)
        end
        mv(stage, output)
    finally
        isdir(stage) && rm(stage; recursive=true)
    end
    println("Saved CGPM and SPMM examples to $output")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
