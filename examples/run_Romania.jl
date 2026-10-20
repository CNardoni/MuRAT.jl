# Complete MuRAT run on the Romania example (1D LITHO1.0 model, needs NCDatasets).
#
# The data (sac_Romania/, velocity_models/) are in this folder. Start Julia with
# threads from the folder where the results should go, and run:
#
#     julia -t auto --project=/path/to/myenv /path/to/MuRAT.jl/examples/run_Romania.jl
#
# or, from the Julia REPL: include(joinpath(pkgdir(MuRAT), "examples", "run_Romania.jl"))
#
# where /path/to/myenv is a Julia environment with MuRAT and the packages loaded
# below, created once with:
#
#     julia --project=/path/to/myenv -e 'using Pkg; Pkg.add(url = "https://github.com/CNardoni/MuRAT.jl"); Pkg.add(["CairoMakie", "NCDatasets"])'
#
# Results go to ./Romania in that folder (tables in Romania/TXT, ParaView files in Romania/Results,
# figures in the usual MuRAT folders, intermediate results in Romania/temp).

using MuRAT
using NCDatasets          # needed: reads the LITHO1.0 velocity model
using CairoMakie          # optional: remove this line to skip the figures

input = include(joinpath(@__DIR__, "input_Romania.jl"))
murat = run_murat(input)

# Results are in the same places as in the MATLAB structure, e.g.
#   murat.Qc.modvQc              Qc model, checkerboard and spike (cells x 9 x frequencies)
#   murat.Q.modvQ                Q model and tests
#   murat.PD.modvPeakDelay       peak-delay variations and hit counts
#   murat.Qc.averageQcFrequency  average Qc^-1 per frequency
#
# A saved run can be reloaded later with
#   murat = murat_load("Romania/temp/Murat_inverse.jld2")
