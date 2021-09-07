include("Load_LTSpice_Net.jl")
include("SPICE2Matrix.jl")


"""
This function runs an AC analysis (freq. sweep) for an LTSPICE netlist file.

    if no inputs are given, it will prompt for the user to pick a file
    
    Inputs:
    FileName: Path to netlist (LTSPICE: View→SPICE netlist, and it will show a text file as well as populate in PWD)
    FreqList: List of frequencies to test (Hz), Default is 100 Hz to 100kHz in steps of 10
    inputs: vector of independent nodal current inputs (1: NumNodes) and independent V inputs (NumNodes+1:End)
            All voltage inputs default to 1
"""
function RunACAnalysis(FileName=nothing, FreqList = 100:10:100e3, inputs = nothing)
    if FileName===nothing
            FileName = open_dialog("Pick a file")
    end
    
    SPICE_DF,NodeList,InputList,NumVSources,G,BL,Bc,ESR_L,ESR_C,SrcMat = SPICE2Matrix(FileName)
    println(InputList)
    if inputs===nothing
        println("No inputs given")
        inputs = zeros(length(InputList))
        inputs[end-(NumVSources-1):end] .= 1
    end
    
    ResultNodeNames = vcat("V(".*NodeList.*")", "I(".*(InputList[end-(NumVSources-1):end]).*")")
    
    
    Results =[ abs.(inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*FreqList[i],InputList))*inputs) for i in 1:length(FreqList)]
    Results = hcat(Results...)
    ResDict = Dict(ResultNodeNames[1] => Results[1,:])
    for i in 2:length(ResultNodeNames)
      
         merge!(ResDict,Dict(ResultNodeNames[i] => Results[i,:]))
       
    end
    
    


    return FreqList,Results, ResDict

end

function DetermineTempCo(FileName=nothing, DriveFreq = 25e3, ComponentName = "Ldrive";
    FreqList = 100:10:100e3,
    inputs = nothing,
    Θᵣ = 1, Θc = 1, # ᵒC per Watt
    TCᵣ = 0.004, TCc = 0.0003,# Ratio drift (UL) per ᵒC
    InputScaling = 1) 



    if FileName===nothing
        FileName = open_dialog("Pick a file")
    end

    SPICE_DF,NodeList,InputList,NumVSources,G,BL,Bc,ESR_L,ESR_C,SrcMat = SPICE2Matrix(FileName)
    println(InputList)
    if inputs===nothing
        println("No inputs given")
        inputs = zeros(length(InputList))
        inputs[end-(NumVSources-1):end] .= (1*InputScaling)
    end

    ResultNodeNames = vcat("V(".*NodeList.*")", "I(".*(InputList[end-(NumVSources-1):end]).*")")


    Results = [ abs.(inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*FreqList[i],InputList))*inputs) for i in 1:length(FreqList)]
    Results = hcat(Results...)
    ResDict = Dict(ResultNodeNames[1] => Results[1,:])
    for i in 2:length(ResultNodeNames)
    
        merge!(ResDict,Dict(ResultNodeNames[i] => Results[i,:]))
    
    end
    CurVec = plotACElCurrent(SPICE_DF,FreqList,Results,ComponentName,inputs)
    plot(FreqList,abs.(CurVec))
    CurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
    CurrResults = hcat(CurrResults...)

    CurrentDF = getElementCurrents(SPICE_DF,CurrResults,DriveFreq)
    CurrentElIndex = findfirst(isequal(ComponentName),CurrentDF.Name)

    BaselineCurrent = CurrentDF.Current[CurrentElIndex]

   
    SPICE_DF.Drift = zeros(length(SPICE_DF.Name))
    SPICE_DF.TempRise = zeros(length(SPICE_DF.Name))
    SPICE_DF.Dissipation = zeros(length(SPICE_DF.Name))
    for i in 1:length(SPICE_DF.Name)
        
        if SPICE_DF.Type[i] =='R'
            Z = SPICE_DF.Value[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²R
            SPICE_DF.TempRise[i]    = Θᵣ*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ) 
            SPICE_DF.Value[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ) * SPICE_DF.Value[i]

        elseif SPICE_DF.Type[i] =='L'
            Z = SPICE_DF.ESR[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²(ESR)
            SPICE_DF.TempRise[i]    = Θᵣ*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ) 
            SPICE_DF.ESR[i]       = (1 + Θᵣ*SPICE_DF.Dissipation[i]*TCᵣ) * SPICE_DF.ESR[i]# for inductors, only the ESR will drift, not inductance
        elseif SPICE_DF.Type[i] =='C'
            Z = SPICE_DF.ESR[i]
            TmpIndex = findfirst(isequal(SPICE_DF.Name[i]),CurrentDF.Name)

            SPICE_DF.Dissipation[i] = abs(CurrentDF.Current[TmpIndex])^2*Z #Power = I²(ESR)
            SPICE_DF.TempRise[i]    = Θc*SPICE_DF.Dissipation[i]
            SPICE_DF.Drift[i]       = (1 + Θc*SPICE_DF.Dissipation[i]*TCc)
            SPICE_DF.Value[i]       = (1 + Θc*SPICE_DF.Dissipation[i]*TCc) * SPICE_DF.Value[i]

        end
    end
    
    NewCurrResults =[ inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*DriveFreq,InputList))*inputs]
    NewCurrResults = hcat(NewCurrResults...)

    NewCurrentDF = getElementCurrents(SPICE_DF,NewCurrResults,DriveFreq)
    PostHeatCurrent = NewCurrentDF.Current[CurrentElIndex]

    MagDriftPercent = abs((PostHeatCurrent - BaselineCurrent) / PostHeatCurrent)*100

    Results2 =[ abs.(inv(SPICE_DF2Matrix_ω(SPICE_DF,2*π*FreqList[i],InputList))*inputs) for i in 1:length(FreqList)]
    Results2 = hcat(Results2...)
    ResDict2 = Dict(ResultNodeNames[1] => Results2[1,:])
    for i in 2:length(ResultNodeNames)
    
        merge!(ResDict2,Dict(ResultNodeNames[i] => Results2[i,:]))
    
    end
    CurVec2 = plotACElCurrent(SPICE_DF,FreqList,Results2,ComponentName)
    plot(FreqList,abs.(CurVec2))


    return SPICE_DF,MagDriftPercent,NewCurrentDF
end

function getElementCurrents(SPICE_DF,Results,Freq)
    NumPasives = sum(SPICE_DF.Type .!= 'V')
    CurrentDF = DataFrame((Name = Any[], ΔV = Any[], Z = Any[],  Current = Any[]))
    for i in 1:length(SPICE_DF.Name)
        if SPICE_DF.Type[i] !='V'
            N1 = SPICE_DF.Node1[i]
            N2 = SPICE_DF.Node2[i]
            if (N1+N2)<999 #if both are non-GND
                ΔV  = Results[N1] - Results[N2]
            elseif N1==999
                ΔV  = Results[N2]
            else
                ΔV  = Results[N1]
            end
        
        end
        if SPICE_DF.Type[i] =='R'
            Z = SPICE_DF.Value[i]
            Cur = ΔV ./ Z
            push!(CurrentDF,[SPICE_DF.Name[i] ΔV Z Cur])

        elseif SPICE_DF.Type[i] =='L'
            Z = 2 .* π .* Freq.*im.*SPICE_DF.Value[i] .+ SPICE_DF.ESR[i]
            Cur = ΔV ./ Z
            push!(CurrentDF,[SPICE_DF.Name[i] ΔV Z Cur])

        elseif SPICE_DF.Type[i] =='C'
            Z = 1 ./( 2 .* π .* Freq.*im.*SPICE_DF.Value[i]) .+ SPICE_DF.ESR[i]
            Cur = ΔV ./ Z
            push!(CurrentDF,[SPICE_DF.Name[i] ΔV Z Cur])

        end


    end
    return CurrentDF
end

function YMatrix(ω,G,BL,Bc,SrcMat,ESR_L,ESR_C) 
    
    YLMat = complex(zeros(size(G)))
    YCMat = complex(zeros(size(G)))
    
    for i in 1:length(G[:,1]),j in 1:length(G[:,1])
        if BL[i,j]==0
            YL = BL[i,j] # if there is no inductor between nodes, the assumption is it is open, so there is no admittance
        elseif abs(ESR_L[i,j])>0
            YL = BL[i,j] ./ (im.*ω)
            YLMat[i,j] = 1 ./( 1 / YL + 1/ESR_L[i,j])
        else
            YLMat[i,j] = BL[i,j] ./ (im.*ω)
        end

        if Bc[i,j]==0
            YC = Bc[i,j]
        elseif abs(ESR_C[i,j])>0
            YC =  im.*ω.* Bc[i,j]
            YCMat[i,j] = 1 ./( 1 ./ YC + 1/ESR_C[i,j])
        else
            YCMat[i,j] = im.*ω.* Bc[i,j]
        end
        
    end

    return G  .+ YLMat .+ YCMat .+  SrcMat
end
    
function plotACElCurrent(SPICE_DF,FreqList,Results,ComponentName)
    ElIndex = findfirst(isequal(ComponentName),SPICE_DF.Name)


        

            N1 = SPICE_DF.Node1[ElIndex]
            N2 = SPICE_DF.Node2[ElIndex]
            if (N1+N2)<999 #if both are non-GND
                ΔV  = Results[N1,:] - Results[N2,:]
            elseif N1==999
                ΔV  = Results[N2,:]
            else
                ΔV  = Results[N1,:]
            end
        
       
        if SPICE_DF.Type[ElIndex] =='R'
            Z = SPICE_DF.Value[ElIndex]
            CurrentVec = ΔV ./ Z
            

        elseif SPICE_DF.Type[ElIndex] =='L'
            Z = [2 .* π .* FreqList[i].*im.*SPICE_DF.Value[ElIndex] .+ SPICE_DF.ESR[ElIndex] for i in 1:length(FreqList)]
            CurrentVec = ΔV ./ Z
            

        elseif SPICE_DF.Type[ElIndex] =='C'
            Z = [1 ./( 2 .* π .* FreqList[i].*im.*SPICE_DF.Value[ElIndex]) .+ SPICE_DF.ESR[ElIndex] for i in 1:length(FreqList)]
            CurrentVec = ΔV ./ Z
            

        end


    
    
    return CurrentVec
end