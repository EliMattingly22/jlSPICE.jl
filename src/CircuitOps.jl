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



"""

Calculates parallel impedances.
If given two inputs, it assumes both are inputs
If given an array, it takes the parallel combo of it all

Eg.
ZList = [1 + 1*im, 1]
Par(ZList[1], ZList[2])
Par(ZList)
"""
function Par(Z1, Z2)
    return 1 / (1 / Z1 + 1 / Z2)
end
function Par(Z::Array)
    Y = 1 ./Z
    YTotal = sum(Y)
    Zeff = 1/YTotal
    return Zeff
end

"""
Parallel component lists
Input is Array{Tuple{Number,String},1}
e.g.
Freq = 25e3
CompList = [(3, "R"),(30e-6, "L"),(200e-9,"C")]
Par(CompList,Freq)

"""
function Par(C::Array{Tuple{Real,String},1},freq)

    Z = Complex.(zeros(length(C)))
    for i in 1:length(C)
        if (lowercase(C[i][2])=="r") | (lowercase(C[i][2])=="resistor")
            Z[i] = C[i][1]
        elseif (lowercase(C[i][2])=="l") | (lowercase(C[i][2])=="inductor" )
            Z[i] = Z_Ind(C[i][1],freq)
        elseif (lowercase(C[i][2])=="c") | (lowercase(C[i][2])=="capacitor" )
            Z[i] = Z_Cap(C[i][1],freq)
        elseif lowercase(C[i][2])=="z"
            Z[i] = C[i][1]
        else
            error("Unknown component")
        end

    end


    Y = 1 ./Z
    YTotal = sum(Y)
    Zeff = 1/YTotal
    return Zeff
end
