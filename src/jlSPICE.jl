module jlSPICE
using DataFrames, Gtk, PyPlot

include("RunAnalysis.jl")
include("Load_LTSpice_Net.jl")
include("SPICE2Matrix.jl")

export RunACAnalysis,
DetermineTempCo,
SPICE2Matrix,
getElementCurrents,
plotACElCurrent,
ProcessSPICE_DF,
SPICE_DF2Matrix_ω,
LTSpiceLoad,
MakeNumericalVals

# Write your package code here.

end
