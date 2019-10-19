# TODO test jldoctests
function genGauss(aG, mG, AG;
                  xmin=nothing, xmax=nothing, nPts=nothing)
    """
    Calculate sum of Gaussians given amplitudes, means, and covariance matrices.

    Evaluates nGaussian Gaussian functions on a grid, returns their sum.
    The grid models space with nPts discretization points along nDims dimensions.
    LaTeX pseudocode:
    
       sum_i aG_i * exp( - (x-mG_i)' * AG_i * (x-mG_i)).
    

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
    size(aG) == (nGauss,) || throw(DimensionMismatch("size of aG, $(size(aG)), does not equal (size(mG,2),), $((size(mG,2),))"))
    size(AG) == (nDims, nDims, nGauss) || throw(DimensionMismatch("size of AG, $(size(AG)), does not equal (size(mG,1),size(mG,1),size(mG,2)), $((size(mG,1),size(mG,1),(size(mG,2))))"))

    # instantiate keyword args
    xmin == nothing ? xmin = - ones(Float64, nDims) : nothing
    xmax == nothing ? xmax = + ones(Float64, nDims) : nothing
    nPts == nothing ? nPts = 67 .* ones(Int64, nDims) : nothing

    # One-dimensional case
    if nDims == 1
        x1 = range(xmin[1], stop=xmax[1], length=nPts[1]);
        g = zeros(nPts[1]);

        for iG = 1:nGauss
            g = g + aG[iG] .* 
                exp.( - (x1 .- mG[iG]) .^2 * AG[iG] .^2)
        end
    end
    
    # Two-dimensional case
    if nDims == 2
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

    if nDims == 3
        x1 = range(xmin[1], stop=xmax[1], length=nPts[1]);
        x2 = range(xmin[2], stop=xmax[2], length=nPts[2]);
        x3 = range(xmin[3], stop=xmax[3], length=nPts[3]);
        g = zeros(nPts[1], nPts[2], nPts[3]);

        for iG = 1:nGauss
            for ix3 = 1:nPts[3], ix2 = 1:nPts[2], ix1 = 1:nPts[1]
                x = [x1[ix1]; x2[ix2]; x3[ix3]]
                g[ix1,ix2,ix3] = g[ix1,ix2,ix3] +
                    aG[iG] * exp.(- transpose(x .- mG[:,iG]) *
                                  AG[:,:,iG] * (x .- mG[:,iG]))
            end
        end
    end

    return g
end


function genRandGauss(; nGauss=30,
                      nPts = [67; 67],
                      aGmin = 0.2,
                      aGmax = 2.0,
                      mGrad = 0.8,
                      AGmin = 10,
                      AGmax = 200,
                      randSeed = nothing,
                      suppShape = "ball")
    """
    Evaluate sum of Gaussians with random amplitudes, means, and covariance matrices.

    Generates random amplitudes, means, and covariance matrices, and calls
    genGauss function. Dimension is determined from the size of nPts arguments.

    # Keyword Arguments
    - `nGauss::Int64`: number of Gaussians
    - `nPts::Array{Int64}`: number of discr. pts along each dimension
    - `aGmin::Float64`
    - `aGmax::Float64` : 
          amplitudes are uniformely distributed between these values.
    - `mGrad::Float64` : 
          means (each coordinate) are  uniformely distributed between 
          -mGrad and mGrad, if suppShape = "box"
    - `mGrad::Float64` : 
          means are  uniformely distributed in the nDim-dimensional 
          ball with radius mGrad, if suppShape = "ball"
    - `AGmin::Float64`
    - `AGmax::Float64` : 
          eigenvalues of covariance matrices are uniformely distributed 
          between these values. AGmin small --> slowly decaying gaussians
          (try AGmin = 10, AGmax = 15). AGmin large --> fast decaying gaussians
          (try AGmin = 500, AGmax = 1000).
    - `randSeed::Int64` : set random seed.
    - `suppShape::String` : "ball" or "box", see mGrad variable above.

    # Examples
    ```jldoctest
    julia> # Default 2D density with 30 Gaussians
    julia> g = genRandGauss();
    julia> # With specific random seed:
    julia> g = genRandGauss(randSeed = 0);
    julia> # The same with better resolution:
    julia> g = genRandGauss(randSeed = 0, nPts=[512,512]);
    julia> # One-dimensional, 10 Gaussians, smaller support
    julia> FPRS.genRandGauss(nGauss = 10, nPts=[512], AGmin=30, AGmax=45, mGrad=0.5, randSeed=0)
    julia> # 50 Gaussians, three-dimensional, with means in a box and faster decay
    julia> g = FPRS.genRandGauss(nGauss = 50, nPts=[80,80,80], AGmin=100, AGmax=1000, mGrad=0.5, randSeed=0, suppShape="box")
    julia> # E.g., plot after integrating with >>> heatmap(sum(g,dims=3)[:,:,1])
    ```
    """

    if suppShape != "box" && suppShape != "ball"
         throw(ArgumentError("suppShape, $(suppShape), must be either 'box' or 'ball'"))
    end
    
    # Convenience variables
    nDims = length(nPts)

    if randSeed == nothing
        seedVal = Random.seed!()
    else
        seedVal = Int(randSeed)
        Random.seed!(seedVal)
    end

    aG = (aGmax-aGmin) * rand(nGauss) .+ aGmin

    if nDims == 1
        mG = zeros(1,nGauss);
        mG[1,:] = 2.0 * mGrad * rand(nGauss) .- mGrad
        AG = zeros(1,1,nGauss);
        AG[1,1,:] = (AGmax-AGmin) * rand(nGauss) .+ AGmin
    end

    if nDims == 2
        # Generate means of Gaussians
        if suppShape == "box"
            mG = zeros(2,nGauss)
            mG[1,:] = 2.0 * mGrad * rand(nGauss) .- mGrad
            mG[2,:] = 2.0 * mGrad * rand(nGauss) .- mGrad
        elseif suppShape == "ball"
            mG = mGrad * randBall(2,nGauss)
        end
        
        # Generate variance matrices of Gaussians
        AG = zeros(2,2,nGauss);
        for iG = 1:nGauss
            # Prepare random rotation matrix (SO(2))
            theta = 2. * pi * rand();
            Omx = [cos(theta) -sin(theta); sin(theta) cos(theta)]
            Dmx = Diagonal((AGmax-AGmin) * rand(2) .+ AGmin)
            Amx = transpose(Omx) * Dmx * Omx
            AG[:,:,iG] = Amx[:,:]
        end
    end

    if nDims == 3
        # Generate means of Gaussians
        if suppShape == "box"
            mG = zeros(3,nGauss)
            mG[1,:] = 2.0 * mGrad * rand(nGauss) .- mGrad
            mG[2,:] = 2.0 * mGrad * rand(nGauss) .- mGrad
            mG[3,:] = 2.0 * mGrad * rand(nGauss) .- mGrad
        elseif suppShape == "ball"
            mG = mGrad * randBall(3,nGauss)
        end
        
        # Generate variance matrices of Gaussians
        AG = zeros(3,3,nGauss);
        for iG = 1:nGauss
            # Prepare random rotation matrix (SO(2))
            Omx = randRotationMatrix()
            Dmx = Diagonal((AGmax-AGmin) * rand(3) .+ AGmin)
            Amx = transpose(Omx) * Dmx * Omx
            AG[:,:,iG] = Amx[:,:]
        end
    end
    
    return genGauss(aG, mG, AG; nPts=nPts)
end
