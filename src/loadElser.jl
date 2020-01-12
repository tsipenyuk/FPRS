"""
    loadElser(fname)

Load (and symmetrize) benchmark phase retrieval data

# Implementation
Loads data from
https://github.com/veitelser/phase-retrieval-benchmarks
provided file name (+path from cwd)

Symmetrize the data (make it twice is large; ignore optimization
due to real-valuedness).

# Examples
```julia-repl
julia> using FPRS;
julia> sqrtI, nAtoms, suppSize = loadElser("../data/data100E")
```
"""
function loadElser(fname::String)
    fIn = readdlm(fname,'\t')
    nAtoms = parse(Int16, fname[end-3:end-1])
    suppSize = 8 * nAtoms
    M = 128
    sqrtI = zeros(M,M)
    sqrtI[:,1:64] = fIn[:,:]
    sqrtI[:,end:-1:66] = fIn[:,2:end]
    return sqrtI, nAtoms, suppSize
end

function loadElserPure(fname::String)
    fInOriginal = readdlm(fname,'\t')
    fIn = sqrt.(fInOriginal) # to match c++ and julia FFTW normalization constants
    fIn = transpose(fIn) # To match julia's rfft format
    fmag[1:64,:] = fIn
    return fmag
end



