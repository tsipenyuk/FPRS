# randBall
@testset "randBall" begin
           res = randBall(2,10)
           @test typeof(res) <: Array{Float64,2}
end
