"""
    randBall(nDims, nPts)

Draw nPts random points uniformely distributed in an nDims-dimensional ball.

# Implementation
Samples uniformely in a unit cube, then throws away unnecessary points.
Efficient only for low number of dimensions (nDims=1,2,3).

# Examples
```julia-repl
julia> using FPRS; using Plot;
julia> x = randBall(2,50);
julia> figure(); scatter(x[1,:], x[2,:])
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
