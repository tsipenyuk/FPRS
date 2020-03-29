module FPRS

using FFTW
using LinearAlgebra
using Random         # for densities
using DelimitedFiles # readdlm for loadElser

include("utilities.jl")
include("densities.jl")
include("projections.jl")
include("loadElser.jl")
include("bookkeeping.jl")

export
# Utilities
randBall
randRotationMatrix
# Densities
genGaussWGrid
genGauss
genRandGaussWGrid
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
loadElserPure
# bookkeeping
pri
pri1D
lpri

end # module
