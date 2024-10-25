function Noise(Zₗ,en,in;T=300)
    k_b = 1.38e-23
     sqrt.(4*k_b*T*real(Zₗ)) .+ en .+ in*abs(Zₗ)
 end


function ThermalNoise(R_L)
    k_b = 1.38e-23
     sqrt.(4*k_b*300*R_L)
end


""" Enter Wire AWG value, returns Diam, Area,MaxFreq"""
function AWG(AWG_Num)
        
    Diam = .005 * 92^((36-AWG_Num)/39)*25.4 #mm
            
    Area = pi/4*Diam^2      
    μ₀ = 4*π*1e-7
    ρ = 1.7e-8 #Ohm-meter
    MaxFreq =  ρ /((Diam*1e-3/2)^2 *μ₀*π)
    return Diam, Area,MaxFreq
end

""" Enter Wire AWG value and length in meters, returns Diam[mm], Area[mm^2],MaxFreq[Hz], DC resistance[Ohms]"""
function AWG(AWG_Num,L)
        
    Diam = .005 * 92^((36-AWG_Num)/39)*25.4 #mm
            
    Area = pi/4*Diam^2      
    μ₀ = 4*π*1e-7
    ρ = 1.7e-8 #Ohm-meter
    MaxFreq =  ρ /((Diam*1e-3/2)^2 *μ₀*π)
    return Diam, Area,MaxFreq,  resistance(L,Area*1e-6)
end

function resistance(L,A;ρ = 1.7e-8) #Ohm-meter)

    return ρ*L/A
end

