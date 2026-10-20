# Complete MuRAT run on the Toba example (1D model iasp91, 3 components).
#
# The data (sac_Toba/, velocity_models/) are in this folder. Start Julia with
# threads from the folder where the results should go, and run:
#
#     julia -t auto --project=/path/to/myenv /path/to/MuRAT.jl/examples/run_Toba.jl
#
# or, from the Julia REPL: include(joinpath(pkgdir(MuRAT), "examples", "run_Toba.jl"))
#
# where /path/to/myenv is a Julia environment with MuRAT and the packages loaded
# below, created once with:
#
#     julia --project=/path/to/myenv -e 'using Pkg; Pkg.add(url = "https://github.com/CNardoni/MuRAT.jl"); Pkg.add(["CairoMakie"])'
#
# Results go to ./Toba in that folder (tables in Toba/TXT, ParaView files in Toba/Results,
# figures in the usual MuRAT folders, intermediate results in Toba/temp).

using MuRAT
using CairoMakie          # optional: remove this line to skip the figures

input = include(joinpath(@__DIR__, "input_Toba.jl"))
murat = run_murat(input)

# Results are in the same places as in the MATLAB structure, e.g.
#   murat.Qc.modvQc              Qc model, checkerboard and spike (cells x 9 x frequencies)
#   murat.Q.modvQ                Q model and tests
#   murat.PD.modvPeakDelay       peak-delay variations and hit counts
#   murat.Qc.averageQcFrequency  average Qc^-1 per frequency
#
# A saved run can be reloaded later with
#   murat = murat_load("Toba/temp/Murat_inverse.jld2")
