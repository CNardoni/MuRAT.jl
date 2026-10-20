# MuRAT input settings and the Murat result container.

"""
    MuratInput(; kwargs...)

User settings of a MuRAT run. Field names are the same as in the MATLAB input
files (`Murat.input.<name>`); the only exception is `end`, a reserved word in
Julia, which is called `ending`. Empty MATLAB values (`[]`) are `nothing`.
See `examples/input_MSH.jl` for a complete input.
MATLAB: the `Murat.input` structure filled by the input script (e.g. `Murat_inputMSH.m`).
Fields without a default are mandatory keywords.

# Fields (units as in the MATLAB input file)
- `dataDirectory`: folder with the SAC files; `label`: output folder name.
- `workers`: kept for compatibility with the MATLAB input files (MATLAB: size of
  the parallel pool). Julia cannot add threads to a running session: start Julia
  with them, e.g. `julia -t 6`.
- `originTime`, `PTime`, `STime`: SAC header fields with origin, P and S times
  (`"o"`, `"a"`, `"t0"` ... `"t9"`); `POrS`: phase used, 2 = P, 3 = S.
- `centralFrequency`: central frequencies (Hz); `envelopeSmoothTime`: envelope
  smoothing time (s); `components`: 1, 2 or 3 components (order E, N, Z);
  `declustering`: grid-division factor for declustering (`nothing` = no declustering).
- `minimumPeakDelay`, `maximumPeakDelay` (s); `stDevPD`: number of standard
  deviations used to keep peak delays close to the average log trend.
- `bodyWindow`: direct-wave (and noise) window length (s); `startNoise`: start of
  the noise window (s); `tresholdNoise`: minimum coda-to-noise energy ratio.
- `lapseTimeMethod`: `"Constant"`, `"Peak"` or `"Travel"`; `startLapseTime`: coda
  start (s after origin) for `"Constant"`, travel-time multiplier for `"Travel"`;
  `codaWindow`: coda window length (s); `mintravel`, `maxtravel`: accepted travel
  times (s); `albedo`, `iExtinctionLength` (1/km): one value per frequency;
  `kernelTreshold`: factor dividing the grid for the kernel computation;
  `kernelType`: `"differential"` (default, MuRAT 4.0: kernel at the coda start minus
  kernel at the coda end) or `"middle"` (MuRAT 3: kernel at the middle of the coda
  window);
  `QcMeasurement`: `"Linearized"` or `"NonLinear"`; `fitTresholdLinear`: threshold
  on the Qc uncertainty used to retain Qc measurements (`Murat_retainQc`).
- `origin`, `ending`: grid corners `[lat, lon, altitude in m]`; `gridLat`, `gridLong`,
  `gridZ`: number of nodes; `sections`: figure sections `[lat, lon, altitude in m]`.
- `availableVelocity`: 0 = 1D model, 1 = 3D model; `namev`: model file name;
  `averageVelocityP`, `averageVelocityS`: average velocities (km/s).
- `inversionMethod` (`"Tikhonov"`, ...), `MaximumIterations`, `MaximumStallIterations`,
  `PlotInversion`; `dampingValueQc`, `smoothingValueQc`, `dampingValueQ`,
  `smoothingValueQ`: one value per frequency, `nothing` = automatic.
- `sizeCheck`: checkerboard block size in nodes (2 or 4); `spikeLocationOrigin`,
  `spikeLocationEnd` (`[lat, lon, altitude in m]`), `spikeValue`: spike test.
- `Plot*`: 1 to create the corresponding figures, 0 to skip them.
- `plotInterpolation`: the result maps (peak delay, Qc, Q) and the checkerboard
  input/output are drawn at the exact section positions, interpolated on a grid
  this many times finer (default 3; 1 = flat cells on the nearest grid planes, as
  the hit map, kernels and spikes). MATLAB `Murat_plot.m` uses a fixed `divi = 5`.
- `SaveMatFiles`: 1 to write MATLAB files (needs `using MAT`): the maps of
  `Murat_plot.m` in `<label>/matfiles` and a compact `<label>/Murat_summary.mat`
  (see [`murat_matfiles`](@ref)); 0 (default) to skip them.
"""
Base.@kwdef mutable struct MuratInput
    # Data and labels
    dataDirectory::String
    label::String
    workers::Union{Nothing,Int} = nothing        # threads are used if Julia was started with them
    velocityFolder::String = "velocity_models"
    # Header times and phase
    originTime::Union{Nothing,String} = nothing  # nothing -> "o"
    PTime::String = "a"
    STime::Union{Nothing,String} = nothing       # nothing -> "t0"
    POrS::Int = 2
    # Frequencies and components
    centralFrequency::Vector{Float64}
    envelopeSmoothTime::Float64 = 1.0
    components::Int = 1
    declustering::Union{Nothing,Int} = nothing
    # Peak delay
    minimumPeakDelay::Float64
    maximumPeakDelay::Float64
    stDevPD::Float64 = 2.0
    # Coda normalization
    bodyWindow::Float64
    startNoise::Float64
    tresholdNoise::Float64
    # Coda attenuation
    lapseTimeMethod::String
    startLapseTime::Union{Nothing,Float64} = nothing
    codaWindow::Float64
    mintravel::Float64
    maxtravel::Float64
    albedo::Vector{Float64}
    iExtinctionLength::Vector{Float64}
    kernelTreshold::Float64
    # "differential" (MuRAT 4.0): K(tCoda) - K(tCoda + codaWindow);
    # "middle" (MuRAT 3): K(tCoda + codaWindow/2)
    kernelType::String = "differential"
    QcMeasurement::String
    fitTresholdLinear::Float64
    # Grid
    origin::Vector{Float64}
    ending::Vector{Float64}
    gridLat::Int
    gridLong::Int
    gridZ::Int
    sections::Vector{Float64}
    # Velocity model
    availableVelocity::Int
    namev::String
    averageVelocityP::Float64
    averageVelocityS::Float64
    # Inversion
    inversionMethod::String = "Tikhonov"
    MaximumIterations::Int = 2000
    MaximumStallIterations::Int = 100
    PlotInversion::Int = 0
    dampingValueQc::Union{Nothing,Vector{Float64}} = nothing
    smoothingValueQc::Union{Nothing,Vector{Float64}} = nothing
    dampingValueQ::Union{Nothing,Vector{Float64}} = nothing
    smoothingValueQ::Union{Nothing,Vector{Float64}} = nothing
    # Synthetic tests
    sizeCheck::Int = 2
    spikeLocationOrigin::Union{Nothing,Vector{Float64}} = nothing
    spikeLocationEnd::Union{Nothing,Vector{Float64}} = nothing
    spikeValue::Union{Nothing,Float64} = nothing
    # Figures
    PlotRays::Int = 1
    PlotTests::Int = 1
    PlotResults::Int = 1
    PlotCheckers::Int = 1
    PlotSpikes::Int = 1
    PlotParameters::Int = 1
    # Interpolation factor of the result and checkerboard maps (1 = flat cells)
    plotInterpolation::Int = 3
    # MATLAB files: 1 writes <label>/matfiles (as Murat_plot.m) and
    # <label>/Murat_summary.mat; needs `using MAT`
    SaveMatFiles::Int = 0
end

"""
    MuratSection()

A group of results (`rays`, `PD`, `Qc`, `Q` or `grid`), accessed like the
MATLAB structure: `murat.Qc.inverseQc`.
MATLAB: a sub-structure of `Murat` (e.g. `Murat.Qc`). Fields are stored in a
`Dict{Symbol,Any}`, so new fields can be added at any time, as in MATLAB.
"""
struct MuratSection
    d::Dict{Symbol,Any}
end
"""
    MuratSection() -> MuratSection

Empty section, to be filled with `section.name = value`.
"""
MuratSection() = MuratSection(Dict{Symbol,Any}())
# The property overloads below use `getfield(s, :d)` because `s.d` would call the
# overloaded getproperty itself.
"""
    getproperty(s::MuratSection, k::Symbol)

`s.k`: value of field `k` (error if missing), like MATLAB `Murat.Qc.k`.
"""
Base.getproperty(s::MuratSection, k::Symbol) = getfield(s, :d)[k]
"""
    setproperty!(s::MuratSection, k::Symbol, v)

`s.k = v`: creates or replaces field `k`.
"""
Base.setproperty!(s::MuratSection, k::Symbol, v) = (getfield(s, :d)[k] = v)
"""
    propertynames(s::MuratSection)

Names of the stored fields (used by tab completion), like MATLAB `fieldnames`.
"""
Base.propertynames(s::MuratSection) = Tuple(keys(getfield(s, :d)))
"""
    haskey(s::MuratSection, k::Symbol) -> Bool

True if field `k` exists (MATLAB `isfield`).
"""
Base.haskey(s::MuratSection, k::Symbol) = haskey(getfield(s, :d), k)
"""
    show(io::IO, s::MuratSection)

Compact display listing the field names in alphabetical order.
"""
Base.show(io::IO, s::MuratSection) = print(io, "MuratSection(", join(sort(collect(keys(getfield(s, :d)))), ", "), ")")

"""
    Murat

Container of a MuRAT run, mirroring the MATLAB `Murat` structure:
`input` (user settings), `grid` (derived settings: file list, headers,
grids, velocity models), and the results `rays`, `PD`, `Qc` and `Q`.
MATLAB: the `Murat` structure. Difference: the derived settings that MATLAB adds to
`Murat.input` (e.g. `Murat.input.modv`, `listSac`, `x`) are in `grid` here, so that
`input` contains only the user settings.
"""
struct Murat
    input::MuratInput
    grid::MuratSection
    rays::MuratSection
    PD::MuratSection
    Qc::MuratSection
    Q::MuratSection
end
"""
    Murat(input::MuratInput) -> Murat

New run container with the given settings and empty sections.
"""
Murat(input::MuratInput) = Murat(input, MuratSection(), MuratSection(), MuratSection(),
                                 MuratSection(), MuratSection())

"""
    outdir(m::Murat, parts...) -> String

Path inside the output folder of a run: `joinpath(m.input.label, parts...)`, e.g.
`outdir(m, "TXT", "file.txt")`. No MATLAB counterpart (Julia helper).
"""
outdir(m::Murat, parts...) = joinpath(m.input.label, parts...)
