"""
Structures to keep information sorted; for 2D arrays
"""
mutable struct pri
    g::Array{Float64,2}     # approximation 
    gInit::Array{Float64,2} # initial state
    # Energies
    Ep::Array{Float64,1}
    Em::Array{Float64,1}
    E::Array{Float64,1}
    # Additional/optional quantities for debugging etc
    p::Array{Float64,2}
    q::Array{Float64,2}
    d::Array{Float64,1}
    d2::Array{Float64,1}
    name::String
end

pri(gInit, nSteps, name="NoName") =
    pri(gInit, gInit,
        zeros(nSteps), zeros(nSteps), zeros(nSteps),
        gInit, gInit,
        zeros(nSteps), zeros(nSteps),
        name)


mutable struct lpri
    g::Array{Float64,2}     # approximation 
    gInit::Array{Float64,2} # initial state
    # Energies
    Ep::Array{Float64,1}
    Em::Array{Float64,1}
    E::Array{Float64,1}
    # Additional/optional quantities for debugging etc
    p::Array{Float64,2}
    q::Array{Float64,2}
    go::Array{Float64,2}
    qo::Array{Float64,2}
    po::Array{Float64,2}
    d::Array{Float64,1}
    d2::Array{Float64,1}
    d3::Array{Float64,1}
    d4::Array{Float64,1}
    d5::Array{Float64,1}
    name::String
end

lpri(gInit, nSteps, name="NoName") =
    pri(gInit, gInit,
        zeros(nSteps), zeros(nSteps), zeros(nSteps),
        gInit, gInit, gInit, gInit, gInit,
        zeros(nSteps), zeros(nSteps), zeros(nSteps), zeros(nSteps), zeros(nSteps),
        name)


"""
    Structures to keep information sorted; for 1D arrays
"""
mutable struct pri1D
    g::Array{Float64,1}     # approximation 
    gInit::Array{Float64,1} # initial state
    # Energies
    Ep::Array{Float64,1}
    Em::Array{Float64,1}
    E::Array{Float64,1}
    # Additional/optional quantities for debugging etc
    p::Array{Float64,1}
    q::Array{Float64,1}
    d::Array{Float64,1}
    d2::Array{Float64,1}
    name::String
end

pri1D(gInit, nSteps, name="NoName") =
    pri(gInit, gInit,
        zeros(nSteps), zeros(nSteps), zeros(nSteps),
        gInit, gInit,
        zeros(nSteps), zeros(nSteps),
        name)
