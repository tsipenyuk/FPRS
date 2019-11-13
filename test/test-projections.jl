using FFTW

# projections and energies
@testset "projections" begin
    a = ones(2,2)
    a[1,1] = -1
    a[2,2] = 2
    sqrtI = abs.(fft(a))
    b = zeros(2,2)
    b[1,1] = 1
    
    # projections
    @test FPRS.pM(b,sqrtI) == [2.5 0.5; 0.5 -0.5]
    @test FPRS.pM0(b,sqrtI) == [2.0 0.0; 0.0 -1.0]
    @test FPRS.pP(a) == [0.0 1.0; 1.0 2.0]
    @test FPRS.pTa(a, 1.6) == [0.0 0.0; 0.0 2.0]
    @test FPRS.pTn(a, 1) == [0.0 0.0; 0.0 2.0]
    @test FPRS.pTn(a, 3) == [0.0 1.0; 1.0 2.0]

    # energies
    @test FPRS.eM(b,sqrtI) == 1.0
    @test FPRS.eM0(b,sqrtI) == 1.0
    @test FPRS.eP(a) == 0.5
    @test FPRS.eTa(a, 1.6) == 1.5
    @test FPRS.eTn(a, 1) == 1.5
    @test FPRS.eTn(a, 3) == 0.5
end
