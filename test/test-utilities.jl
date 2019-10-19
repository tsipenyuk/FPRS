# randBall
@testset "randBall" begin
    res = FPRS.randBall(2,10)
    @test typeof(res) <: Array{Float64,2}
end

# randRotationMatrix(deflection=1.0)
@testset "randRotationMatrix" begin
           res = FPRS.randRotationMatrix()
           @test typeof(res) <: Array{Float64,2}
end
