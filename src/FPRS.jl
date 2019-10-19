module FPRS

using FFTW
using LinearAlgebra
using Random

include("utilities.jl")
include("densities.jl")

export
    # Utilities
    randBall
    randRotationMatrix
    # Densities
    genGauss

end # module
