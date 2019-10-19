# genGauss
@testset "genGauss" begin
    # One-dimensional test
    aG = ones(3);
    mG = zeros(1,3);
    mG[1,:] = [0.0, 0.7, -0.3];
    AG = zeros(1,1,3);
    AG[1,1,:] = [100 80 120];
    g = FPRS.genGauss(aG, mG, AG);
    @test abs.(sum(g) - 2.8286690853954575) < 67.0 * eps()
    
    # Two-dimensional test
    # Generate 3 Gaussians
    aG = ones(3);
    mG = zeros(2,3);
    mG[1,:] = [0.0, 0.7, -0.3];
    mG[2,:] = [0.0, -0.6, -0.06];
    AG = zeros(2,2,3);
    AG[:,:,1] = [100 0.05; 0 100];
    AG[:,:,2] = [100 0.15; 0 100];
    AG[:,:,3] = [100 -0.3; 0.4 100];
    g = FPRS.genGauss(aG, mG, AG);
    @test abs.(sum(g) - 102.63572536664947) < 67.0^2 * eps()

    xmin = [-1.0 -2.0];
    xmax = [2.0   1.0];
    nPts = [127 127];
    g = FPRS.genGauss(aG, mG, AG, nPts = nPts, xmin=xmin, xmax=xmax);
    @test abs.(sum(g) - 166.25310747322118) < 127.0^2 * eps()
    
    # Test errors
    aG = ones(2)
    @test_throws DimensionMismatch FPRS.genGauss(aG, mG, AG);
    aG = ones(3)
    AG = zeros(2,3,3)
    @test_throws DimensionMismatch FPRS.genGauss(aG, mG, AG);

    # Three-dimensional test
    aG = ones(3);
    mG = zeros(3,3);
    mG[1,:] = [0.0, 0.7, -0.3];
    mG[2,:] = [0.0, -0.6, -0.06];
    mG[3,:] = [0.3, -0.5, -0.2];
    AG = zeros(3,3,3);
    AG[:,:,1] = [100 0.05 50; 0   100 20; 0  0  100];
    AG[:,:,2] = [100 0.15 10; 0   100 10; 0  20 100];
    AG[:,:,3] = [100 -0.3 10; 0.4 100 50; 50 0  100];
    g = FPRS.genGauss(aG, mG, AG);
    @test abs.(sum(g) - 627.7968551738583) < 67.0^3 * eps()
end

# genRandGauss
@testset "genRandGauss" begin
    # 1D
    g = FPRS.genRandGauss(nGauss = 50, nPts=[512], AGmin=100, AGmax=1000, mGrad=0.5, randSeed=0, suppShape="box")
    @test abs.(sum(g) - 63.77852000489939) < 512 * eps()
    
    # 2D
    g = FPRS.genRandGauss(nGauss = 10, nPts=[51,51], mGrad=0.5, randSeed=0, suppShape="ball")
    @test abs.(sum(g) - 210.31303142418054) < 51^2 * eps()

    # 3D
    g = FPRS.genRandGauss(nGauss = 10, nPts=[31,31,21], randSeed=0, suppShape="ball")
    @test abs.(sum(g) - 200.71210294764373) < 31^2 * 21 * eps()

    # Reverse-engineered. (?) Why not ArgumentError 
    @test_throws MethodError FPRS.genGauss(suppShape = 1);
end
