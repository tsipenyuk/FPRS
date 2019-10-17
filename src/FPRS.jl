module FPRS

using FFTW
using LinearAlgebra
using Random

export
    # Utilities
    randBall

include("utilities.jl")
include("densities.jl")

end # module
