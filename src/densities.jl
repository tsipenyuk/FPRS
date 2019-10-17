


# nDims
#  Integer
#  number of dimensions (1,2,3)
#
# nPts
#  vector of integers
#  number of discretization points (resolution) along each dimension
#
# nGauss
#  integer
#  number of Gaussians
#
# mG[:, iGauss]
#  vector of size nDims
#  mean of i-th Gaussian
#
# AG[:, :, iGauss]
#  array of size nDims x nDims
#  inverse of the covariance matrix of i-th Gaussian
#
# xmin xmax
# xG aG vG mG
function genGauss(aG, mG, AG;
                  xmin=nothing, xmax=nothing, nPts=nothing)

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
        x2 = range(xmax[2], stop=xmax[2], length=nPts[2]);
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

nGauss = 3
xmin = -ones(2)
xmax = ones(2)
nPts = 67 * ones(Int64, 2)
aG = ones(3)
mG = zeros(2,3)
x1 = range(xmin[1], stop=xmax[1], length=nPts[1]);
x2 = range(xmax[2], stop=xmax[2], length=nPts[2]);
g = zeros(nPts);

mG[1,:] = [0.1, 0.7, -0.3]
mG[2,:] = [-0.4, -0.6, -0.06]

AG = zeros(2,2,3)
AG[:,:,1] = [1 0.5; 0 1]
AG[:,:,2] = [1 1.5; 0 1]
AG[:,:,3] = [1 -0.5; 2.5 1]

g = zeros(nPts[1],nPts[2]);
for iG = 1:nGauss
    for ix2 = 1:nPts[2], ix1 = 1:nPts[1]
        x = [x1[ix1]; x2[ix2]]
        g[ix1,ix2] = g[ix1,ix2] +
            exp.(- transpose(x-mG[:,iG]) * AG[:,:,iG] * (x-mG[:,iG]))
    end
end

xt = zeros(2, nPts[1], nPts[2])
xt[1,:,:] = x1[:,:]
xt[2,:,:] = x2[:,:]
gt = zeros(nPts[1],nPts[2]);
for iG = 1:nGauss
        x = [x1[ix1]; x2[ix2]]
        g[ix1,ix2] = g[ix1,ix2] +
            exp.(- transpose(x-mG[:,iG]) * AG[:,:,iG] * (x-mG[:,iG]))
end

function test(;x)
    if @isdefined x
        return 1
    else
        return 0
    end
end
