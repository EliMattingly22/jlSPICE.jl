
function BruteForceOpt_MaxFreq(Proc_DF::ProcessedDF,FreqList,Inputs)
# inputs should be an Nx1 array of pairs where first item in pair is the name of the variable, and the second is a 1xM vector of possible values

    VarLengths = zeros(Int8,1,length(Inputs))
    VarIndex = zeros(Int8,size(VarLengths))
    for i in 1:length(Inputs)
        VarLengths[i]=length(Inputs[i][2])
        VarIndex[i]=findfirst(x -> x==Inputs[i][1],convert(Array{String,1},Proc_DF.DF.Name) )

    end
    #FinalData = [[zeros(size(VarLengths))] [zeros(length(Proc_DF.InputList)-1,length(FreqList))]]

    for ii in 1:VarLengths[1], j in 1:VarLengths[2],k in 1:VarLengths[3]

        Proc_DF.DF.Value[VarIndex[1]]=Inputs[1][2][ii]
        Proc_DF.DF.Value[VarIndex[2]]=Inputs[2][2][j]
        Proc_DF.DF.Value[VarIndex[3]]=Inputs[3][2][k]
        TrialRes = RunDF(Proc_DF,FreqList)
        if (ii==1)&(j==1)&(k==1)
            FinalData=TrialRes
            println(FinalData)
        else
            FinalData = [FinalData;TrialRes]
        end
        #FinalData = vcat(FinalData,TrialRes[3,:])

    end
    #plot(log10.(abs.(TrialRes[3,:])))


return FinalData
end

function RunDF(DF::ProcessedDF,FreqList)


    CM = SPICE_DF2Matrix(DF)
    CondMat(w) =
        CM.RMat .+ CM.invLMat ./ (im .* w) .+ im .* w .* CM.CMat .+ CM.SrcMat
    inputs = [0; 0; 0; 0; 0; 1]

    Results = zeros(length(inputs), length(FreqList))
    for i = 1:length(FreqList)
        Freq = FreqList[i]
        Results[:, i] = abs.(inv(CondMat(2 * pi * Freq)) * inputs)
    end
    return Results
end
