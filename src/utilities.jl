"""
    randBall(nDims, nPts)

Draw nPts random points uniformely distributed in an nDims-dimensional ball.

# Implementation
Samples uniformely in a unit cube, then throws away unnecessary points.
Efficient only for low number of dimensions (nDims=1,2,3).

# Examples
```julia-repl
julia> using FPRS; using Plots;
julia> x = randBall(2,50);
julia> figure(); scatter(x[1,:], x[2,:]);
```
"""
function randBall(nDims::Integer, nPts::Integer)

    res = zeros(nDims, nPts)
    numInBall = 0
    
    while numInBall < nPts
        newPoint = 2.0 * rand(nDims) .- 1.0
        if norm(newPoint) <= 1
            numInBall += 1
            res[1:end, numInBall] = newPoint
        end
    end
    
    return res
end


"""
    randRotationMatrix(deflection=1.0; randSeed)

Return a random rotation matrix in 3d.

Sources:
-http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
-http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c

# Arguments
- `deflection`: number between 0 and 1 controls the size of the perturbation. 
                For 0, no rotation; for 1, competely random rotation. 
                Small deflection => small perturbation.

# Examples
```julia-repl
julia> mx = randRotationMatrix();
```
"""
function randRotationMatrix(deflection=1.0)
    # Random.seed!(seedVal) # do not reseed rng
    
    theta = 2.0*deflection* pi * rand() # Rotation about the pole (Z)
    phi = 2.0* pi * rand() # For direction of pole deflection
    z = 2.0*deflection * rand() # For magnitude of pole deflection
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    
    r =  sqrt(z)
    Vx = sin(phi) * r
    Vy = cos(phi) * r
    Vz = sqrt(2.0 - z)

    # Compute the row vector S = Transpose(V) * R, where R is a simple
    # rotation by theta about the z-axis.  No need to compute Sz since
    # it's just Vz.

    st = sin(theta)
    ct = cos(theta)
    Sx = Vx * ct - Vy * st
    Sy = Vx * st + Vy * ct
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R, which
    # is equivalent to V S - R.

    M = zeros(3,3)
    M[1,1] = Vx * Sx - ct;
    M[1,2] = Vx * Sy - st;
    M[1,3] = Vx * Vz;

    M[2,1] = Vy * Sx + st;
    M[2,2] = Vy * Sy - ct;
    M[2,3] = Vy * Vz;

    M[3,1] = Vz * Sx;
    M[3,2] = Vz * Sy;
    M[3,3] = 1.0 - z;   # This equals Vz * Vz - 1.0

    return M
end
