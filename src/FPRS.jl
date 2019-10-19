module FPRS

using FFTW
using LinearAlgebra
using Random

include("utilities.jl")
include("densities.jl")
include("projections.jl")

export
# Utilities
randBall
randRotationMatrix
# Densities
genGauss
genRandGauss
# Projections
pM
pM0
pP
pTa
pTn
eM
eM0
eP
eTa
eTn

end # module
