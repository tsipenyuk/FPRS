module FPRS

using FFTW
using LinearAlgebra
using Random

include("utilities.jl")
include("densities.jl")

export
    # Utilities
    randBall
    # Densities
    genGauss

end # module
