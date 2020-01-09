module FPRS

using FFTW
using LinearAlgebra
using Random         # for densities
using DelimitedFiles # readdlm for loadElser

include("utilities.jl")
include("densities.jl")
include("projections.jl")
include("loadElser.jl")

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
# loadElser
loadElser

end # module
