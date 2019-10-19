"""
Provides a collection of projection operators and error functionals.
"""

###### Projections ######

"Fourier modulus projection"
function pM(g, sqrtI)
    ang = angle.(fft(g))
    res = real(ifft(sqrtI .* exp.(im .* ang)))
end


"""
Incomplete (Elser) modulus projection 

A variant of the modulus projection that leaves 0'th frequency
invariant or sets it to 0 if it is negative. This projection is
used, e.g., in 
> "Benchmark problems for phase retrieval", V. Elser, T.-Y. Lan & T. Bendory
"""
function pM0(g, sqrtI)
    goalAbs = sqrtI
    goalAbs[1,1] = maximum([sum(g), 0.0])
    return real(ifft(goalAbs .* exp.(im .* angle.(fft(g)))))
end


"Positivity projection"
function pP(g)
    res = real(g .* (g .>= 0))
end

"Soft sparsity projection"
function pTa(g, alpha=0)
    res = real(g .* (g .>= alpha))
end


"""
Hard sparsity projection

Set all but suppSize largest elements of g to zero.
"""
function pTn(g, suppSize)
    whereLargest = partialsortperm(g[:], 1:suppSize, rev=true)
    res = zeros(size(g))
    res[whereLargest] = g[:][whereLargest]
    return res
end


###### Error functionals ######
"Modulus energy functional"
function eM(g, sqrtI)
    res = 0.5 * sum((abs.(fft(g)) - sqrtI).^2) / length(g)
end

"Incomplete modulus energy functional"
function eM0(g, sqrtI)
        res = 0.5 * sum((abs.(fft(g)[2:end]) - sqrtI[2:end]).^2) / length(g)
end

"Positive energy functional"
function eP(g)
    res = 0.5 * sum(real(g .* (g .< 0)).^2)
end

"Soft sparsity energy functional"
function eTa(g, alpha=0)
    res = 0.5 * sum((g - pTa(g, alpha)).^2)
end

"Hard sparsity energy functional"
function eTn(g, suppSize)
    res = 0.5 * sum((g - pS0(g, suppSize)).^2)
end

