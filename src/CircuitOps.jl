function Z_Cap(C, f)
    return (1 ./ (im * 2 * pi .* f .* C))
end

function Z_Ind(L, f)
    return im * 2 * pi .* f .* L
end

"""
This function takes in some reactance element (E) and a frequency, and finds the pair


"""
function findResPair(E, f)
    ##Finds the resonant pair for a reactive element, E
    ## ωL = 1/(ωC) → L = 1/(ω^2C),and vice versa
    return 1 ./ ((2 * pi * f) .^ 2 * E)
end