# loadElser
@testset "loadElser" begin
    sqrtI, nAtoms, suppSize = FPRS.loadElser("./testData/testData123X")
    @test typeof(sqrtI) <: Array{Float64,2}
    @test nAtoms == 123
end
