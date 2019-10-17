# TODO write tests
# TODO test jldoctests
function genGauss(aG, mG, AG;
                  xmin=nothing, xmax=nothing, nPts=nothing)
    """
    Calculate sum of Gaussians given amplitudes, means, and covariance matrices.

    Evaluates nGaussian Gaussian functions on a grid, returns their sum. 
    The grid models space with nPts discretization points along nDims dimensions.

    Arguments nGauss and nDims are not declared explicitely, but calculated from
    means vector mG, see below.

    # Arguments
    - `aG::Array{Float64}`: amplitudes of Gaussians, vector of length nGauss
    - `mG::Array{Float64}`: coordinates of Gaussians, 
                            matrix of size nDims x nGauss
    - `AG::Array{Float64}`: covariance matrices of Gaussians, 
                            tensor of size nDims x nDims x nGauss

    # Keyword Arguments
    - `xmin[=-ones(nDims)]`: lower bound of the discretization box
    - `xmax[= ones(nDims)]`: upper bound of the discretization box
    - `nPts[= 67 .* ones(Int64, nDims)]`: number of discretization points
                                          along each dimension

    # Examples
    ```jldoctest
    julia> # Generate 3 Gaussians
    julia> aG = ones(3);
    julia> mG = zeros(2,3);
    julia> mG[1,:] = [0.0, 0.7, -0.3];
    julia> mG[2,:] = [0.0, -0.6, -0.06];
    julia> AG = zeros(2,2,3);
    julia> AG[:,:,1] = [100 0.05; 0 100];
    julia> AG[:,:,2] = [100 0.15; 0 100];
    julia> AG[:,:,3] = [100 -0.3; 0.4 100];
    julia>
    julia> g = FPRS.genGauss(aG, mG, AG);
    julia> g[32:35,33:36]
    3×3 Array{Float64,2}:
     0.912746  0.832567  0.632048
     1.00009   0.912317  0.692624
     0.912275  0.832193  0.631775
    julia> 
    julia> # Double the resolution
    julia> # Make grid bounds larger if necessary
    julia> xmin = [-1.0 -2.0];
    julia> xmax = [2.0   1.0];
    julia> g = FPRS.genGauss(aG, mG, AG, nPts = [127,127], xmin=xmin, xmax=xmax);
    julia> g[71:73, 59:61]
    3×3 Array{Float64,2}:
     0.862873  0.892834  0.824811
     0.955636  0.988733  0.913326
     0.944926  0.977569  0.902937
    julia>
    julia> # At which x1 coordinate was this Gaussian? Calculate with
    julia> x1 = range(xmin[1], stop=xmax[1], length=127);
    julia> collect(x1[73:75])
    3-element Array{Float64,1}:
     0.7142857142857143
     0.7380952380952381
     0.7619047619047619
    ```
    """
    # Convenience variables
    nDims  = size(mG,1)
    nGauss = size(mG,2)

    # Check input size
    size(aG) == (nGauss,) || DimensionMismatch("size of aG, $(size(aG)), does not equal (size(mG,2),), $((size(mG,2),))")
    size(AG) == (nDims, nDims, nGauss) || DimensionMismatch("size of AG, $(size(AG)), does not equal (size(mG,1),size(mG,1),size(mG,2)), $((size(mG,1),size(mG,1),(size(mG,2))))")

    # instantiate keyword args
    xmin == nothing ? xmin = - ones(Float64, nDims) : nothing
    xmax == nothing ? xmax = + ones(Float64, nDims) : nothing
    nPts == nothing ? nPts = 67 .* ones(Int64, nDims) : nothing
        
    # Two-dimensional case
    if length(xmin) == 2
        x1 = range(xmin[1], stop=xmax[1], length=nPts[1]);
        x2 = range(xmin[2], stop=xmax[2], length=nPts[2]);
        g = zeros(nPts[1], nPts[2]);

        for iG = 1:nGauss
            for ix2 = 1:nPts[2], ix1 = 1:nPts[1]
                x = [x1[ix1]; x2[ix2]]
                g[ix1,ix2] = g[ix1,ix2] +
                    aG[iG] * exp.(- transpose(x .- mG[:,iG]) *
                                  AG[:,:,iG] * (x .- mG[:,iG]))
            end
        end
    end
    
    return g
end

#nGauss = 3
#xmin = -ones(2)
#xmax = ones(2)
#nPts = 67 * ones(Int64, 2)
#aG = ones(3)
#mG = zeros(2,3)
#x1 = range(xmin[1], stop=xmax[1], length=nPts[1]);
#x2 = range(xmax[2], stop=xmax[2], length=nPts[2]);
#g = zeros(nPts);
#
#mG[1,:] = [0.1, 0.7, -0.3]
#mG[2,:] = [-0.4, -0.6, -0.06]
#
#AG = zeros(2,2,3)
#AG[:,:,1] = [1 0.5; 0 1]
#AG[:,:,2] = [1 1.5; 0 1]
#AG[:,:,3] = [1 -0.5; 2.5 1]
#
#g = zeros(nPts[1],nPts[2]);
#for iG = 1:nGauss
#    for ix2 = 1:nPts[2], ix1 = 1:nPts[1]
#        x = [x1[ix1]; x2[ix2]]
#        g[ix1,ix2] = g[ix1,ix2] +
#            exp.(- transpose(x-mG[:,iG]) * AG[:,:,iG] * (x-mG[:,iG]))
#    end
#end
#
#xt = zeros(2, nPts[1], nPts[2])
#xt[1,:,:] = x1[:,:]
#xt[2,:,:] = x2[:,:]
#gt = zeros(nPts[1],nPts[2]);
#for iG = 1:nGauss
#        x = [x1[ix1]; x2[ix2]]
#        g[ix1,ix2] = g[ix1,ix2] +
#            exp.(- transpose(x-mG[:,iG]) * AG[:,:,iG] * (x-mG[:,iG]))
#end
#
#function test(;x)
#    if @isdefined x
#        return 1
#    else
#        return 0
#    end
#end
#
