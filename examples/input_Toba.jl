# Toba example (translated from Murat_inputToba.m of MuRAT 4.0)
# Same field names as Murat.input.<name> in MATLAB
# (`end` is called `ending`; MATLAB [] is `nothing`).
# The data are in this folder (examples/); the results go to the folder from
# which the run is started.

using MuRAT

MuratInput(
    # Data folder and output label
    dataDirectory          = joinpath(@__DIR__, "sac_Toba"),
    label                  = "Toba",
    velocityFolder         = joinpath(@__DIR__, "velocity_models"),

    # SAC header times: origin (nothing = "o"), P and S picks (nothing = "t0")
    originTime             = "o",
    PTime                  = "a",
    STime                  = "t0",

    # Phase used for picking: 2 = P, 3 = S
    POrS                   = 3,

    # Central frequencies (Hz) and number of components
    centralFrequency       = [1.5, 3],
    envelopeSmoothTime     = 1,
    components             = 3,

    # Declustering factor (nothing = no declustering)
    declustering           = 5,

    # Peak delay: limits (s) and outlier threshold (standard deviations)
    minimumPeakDelay       = 0.1,
    maximumPeakDelay       = 10,
    stDevPD                = 2,

    # Coda normalization: direct-wave window, noise start (s), noise threshold
    bodyWindow             = 1,
    startNoise             = 3,
    tresholdNoise          = 3,

    # Coda: lapse-time method ("Constant", "Peak" or "Travel"), start, window
    lapseTimeMethod        = "Constant",
    startLapseTime         = 60,
    codaWindow             = 20,
    mintravel              = 12,
    maxtravel              = 35,

    # Kernels: albedo and inverse extinction length per frequency
    albedo                 = [0.5, 0.5],
    iExtinctionLength      = [0.02, 0.02],
    kernelTreshold         = 1,
    # "differential" (MuRAT 4.0, default) or "middle" (MuRAT 3)
    kernelType             = "differential",

    # Qc measurement ("Linearized" or "NonLinear") and fit threshold
    QcMeasurement          = "NonLinear",
    fitTresholdLinear      = 0,

    # Grid: origin and end [lat, lon, altitude (m)], number of nodes
    origin                 = [0, 96, 2000],
    ending                 = [4, 100, -30000],
    gridLat                = 17,
    gridLong               = 17,
    gridZ                  = 11,
    sections               = [2.5, 99, -10000],

    # Velocity model: 0 = 1D, 1 = 3D (file in velocityFolder)
    availableVelocity      = 0,
    namev                  = "iasp91.txt",
    averageVelocityP       = 7,
    averageVelocityS       = 4,

    # Inversion
    inversionMethod        = "Tikhonov",
    MaximumIterations      = 2000,
    MaximumStallIterations = 100,
    PlotInversion          = 0,
    dampingValueQc         = [10, 10],
    smoothingValueQc       = nothing,
    dampingValueQ          = nothing,
    smoothingValueQ        = nothing,

    # Synthetic tests: checkerboard block size and spike [lat, lon, altitude]
    sizeCheck              = 2,
    spikeLocationOrigin    = [2, 98.5, -10000],
    spikeLocationEnd       = [4, 99.5, -20000],
    spikeValue             = 0.02,

    # Figures (1 = yes)
    PlotRays               = 1,
    PlotTests              = 1,
    PlotResults            = 1,
    PlotCheckers           = 1,
    PlotSpikes             = 1,
    PlotParameters         = 1,
)
