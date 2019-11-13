using FPRS
using Test

@testset "FPRS.jl" begin
    include("test-utilities.jl")
    include("test-densities.jl")
    include("test-projections.jl")
end
