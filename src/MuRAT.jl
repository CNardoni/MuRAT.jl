"""
    MuRAT

Julia translation of MuRAT (Multi-Resolution seismic Attenuation
Tomography, https://github.com/LucaDeSiena/MuRAT): peak-delay, coda
attenuation (Qc) and coda-normalization (Q) tomography from SAC data.

Main entry point: [`run_murat`](@ref). The stages can also be run one by
one: [`murat_checks`](@ref), [`murat_data`](@ref) (includes the selection),
[`murat_inversion`](@ref), [`murat_summary`](@ref), [`murat_paraview`](@ref)
and, after `using CairoMakie`, `murat_plot`.

Settings are given as a [`MuratInput`](@ref) (MATLAB `Murat.input`); all results are
stored in a [`Murat`](@ref) container (MATLAB `Murat` structure).
"""
module MuRAT

using DelimitedFiles
using DSP
using FFTW
using JLD2
using LinearAlgebra
using MiniQhull
using Printf
using Random
using SparseArrays
using Statistics
using WriteVTK

export MuratInput, Murat, run_murat
export murat_checks, murat_data, murat_selection, murat_inversion, murat_summary
export murat_paraview, murat_save, murat_save_inversion, murat_load, murat_plot, murat_matfiles
export read_sac, murat_envelope, murat_kernels, murat_coda_matrix, murat_rays

include("utils.jl")
include("geodesy.jl")
include("interpolation.jl")
include("io_sac.jl")
include("input.jl")
include("checks.jl")
include("measurements.jl")
include("rays.jl")
include("kernels.jl")
include("data.jl")
include("selection.jl")
include("inversion.jl")
include("output.jl")

"""
    murat_plot(m::Murat)

Figures of a run (maps, tests, checkerboards, spikes, rays, kernels,
L-curves). Available after `using CairoMakie`.
MATLAB: `Murat_plot.m`. Only the empty generic function is declared here; the method is
defined in the package extension `ext/MuratMakieExt.jl`, loaded with CairoMakie.

# Arguments
- `m`: a run after [`murat_inversion`](@ref) (and [`murat_summary`](@ref)).
"""
function murat_plot end

"""
    murat_matfiles(m::Murat)

MATLAB `.mat` files of a run: the per-frequency maps written by `Murat_plot.m` in
`<label>/matfiles` (peak delay, Qc and Q on the inversion grid and interpolated,
peak-delay hit map, Qc checkerboard input/output; same file and variable names)
and `<label>/Murat_summary.mat`, the `Murat` structure without the largest arrays
(inversion matrices, ray paths, ...), which stays small enough to open in MATLAB.
Available after `using MAT`; called by [`run_murat`](@ref) when the input has
`SaveMatFiles = 1`.
MATLAB: the `save` calls of `Murat_plot.m`; `Murat_summary.mat` replaces the full
`Murat.mat` of `MuRAT.m`, which is several GB for large data sets.
The method is defined in the package extension `ext/MuratMATExt.jl`.

# Arguments
- `m`: a run after [`murat_inversion`](@ref) and [`murat_summary`](@ref).
"""
function murat_matfiles end

"""
    read_litho1(path, lat, lon) -> Matrix

1D velocity profile of the LITHO1.0 model (NetCDF file `path`) at the model node
nearest to (`lat`, `lon`), as a table `[depth (km), radius (km), vp, vs (km/s)]`
every 1 km from 5 km above sea level, used like a 1D model file
(`availableVelocity = 0`, `namev = "LITHO1.0.nc"`).
MATLAB: `Murat_readLithos1.m` (called by `Murat_checks.m` at the grid centre).
The method is defined in the package extension `ext/MuratNCExt.jl`, loaded with
`using NCDatasets`.
"""
function read_litho1 end

"""
    run_murat(input::MuratInput; plots = true, topo = nothing) -> Murat

Complete MuRAT run, as `MuRAT.m`: checks, data processing and selection,
inversions, summary values, tables, ParaView files and (if CairoMakie is
loaded and `plots = true`) figures. Intermediate results are saved in
`<label>/temp` (compressed JLD2): `Murat_checks.jld2`, `Murat_forward.jld2` (the run
after the data stage) and `Murat_inverse.jld2` (only the inversion results; loading
it with [`murat_load`](@ref) adds them to the forward run).
MATLAB: main script `MuRAT.m` (calls `Murat_checks`, `Murat_dataParallelized`,
`Murat_inversion`, `Murat_plot`). Differences: the `temp` files are JLD2, not `.mat`,
are written inside `<label>/temp` instead of `./temp`, and the inverse file does not
repeat the forward data.

# Arguments
- `input`: user settings (MATLAB `Murat.input`).
- `plots`: create the figures (needs `using CairoMakie` beforehand).
- `topo`: optional topography for the vertical sections of the figures, e.g. the
  `GeoData` returned by GeophysicalModelGenerator `import_topo` (see `murat_plot`).

# Returns
- the `Murat` container with all settings and results (also written to
  `<label>/temp/Murat_inverse.jld2`).
"""
function run_murat(input::MuratInput; plots::Bool = true, topo = nothing)
    # Stage 1: folders, SAC header checks, grids and velocity models
    t0 = time()
    m = murat_checks(input)
    murat_save(outdir(m, "temp", "Murat_checks.jld2"), m)
    @printf("Time spent on checks: %.2f s\n", time() - t0)

    # Stage 2: rays, kernels, PD/Qc/Q measurements and data selection
    t1 = time()
    murat_data(m)
    murat_save(outdir(m, "temp", "Murat_forward.jld2"), m)
    @printf("Time spent on data processing: %.2f s\n", time() - t1)

    # Stage 3: inversions and synthetic tests, summary values, ParaView output
    t2 = time()
    murat_inversion(m)
    murat_summary(m)
    # Only the inversion results: the forward data are in Murat_forward.jld2
    murat_save_inversion(outdir(m, "temp", "Murat_inverse.jld2"), m)
    murat_paraview(m)
    @printf("Inversion completed (%.2f s)\n", time() - t2)

    # MATLAB files; the murat_matfiles method exists only once MAT is loaded
    if input.SaveMatFiles == 1
        if hasmethod(murat_matfiles, Tuple{Murat})
            t3 = time()
            murat_matfiles(m)
            @printf("MATLAB files written (%.2f s)\n", time() - t3)
        else
            @info "MATLAB files skipped: run `using MAT` before run_murat to write them."
        end
    end

    # The murat_plot method exists only once the CairoMakie extension is loaded
    if plots
        if hasmethod(murat_plot, Tuple{Murat})
            t4 = time()
            murat_plot(m; topo = topo)
            @printf("Figures completed (%.2f s)\n", time() - t4)
        else
            @info "Figures skipped: run `using CairoMakie` before run_murat to create them."
        end
    end
    @printf("Total time: %.2f s\n", time() - t0)
    return m
end

end # module
