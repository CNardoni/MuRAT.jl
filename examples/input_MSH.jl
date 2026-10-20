# MuRAT input for the Mount St. Helens example, translated from
# Murat_inputMSH.m. Same field names as Murat.input.<name> in MATLAB
# (`end` is called `ending`; MATLAB [] is `nothing`).
# The data are in this folder (examples/); the results go to the folder from
# which the run is started.

using MuRAT

MuratInput(
    # Data folder and output label
    dataDirectory          = joinpath(@__DIR__, "sac_MSH"),
    label                  = "MSH",
    velocityFolder         = joinpath(@__DIR__, "velocity_models"),

    # SAC header times: origin (nothing = "o"), P and S picks (nothing = "t0")
    originTime             = nothing,
    PTime                  = "a",
    STime                  = nothing,

    # Phase used for picking: 2 = P, 3 = S
    POrS                   = 2,

    # Central frequencies (Hz) and number of components
    centralFrequency       = [3, 6, 12, 18],
    envelopeSmoothTime     = 1,
    components             = 1,

    # Declustering factor (nothing = no declustering)
    declustering           = 5,

    # Peak delay: limits (s) and outlier threshold (standard deviations)
    minimumPeakDelay       = 0.1,
    maximumPeakDelay       = 7,
    stDevPD                = 2,

    # Coda normalization: direct-wave window, noise start (s), noise threshold
    bodyWindow             = 1,
    startNoise             = 5,
    tresholdNoise          = 3,

    # Coda: lapse-time method ("Constant", "Peak" or "Travel"), start, window
    lapseTimeMethod        = "Constant",
    startLapseTime         = 15,
    codaWindow             = 10,
    mintravel              = 0.2,
    maxtravel              = 13,

    # Kernels: albedo and inverse extinction length per frequency
    albedo                 = [0.8, 0.6, 0.4, 0.2],
    iExtinctionLength      = [0.12, 0.12, 0.1, 0.1],
    kernelTreshold         = 2,
    # "differential" (MuRAT 4.0, default) or "middle" (MuRAT 3)
    kernelType             = "differential",

    # Qc measurement ("Linearized" or "NonLinear") and fit threshold
    QcMeasurement          = "NonLinear",
    fitTresholdLinear      = 0.006,

    # Grid: origin and end [lat, lon, altitude (m)], number of nodes
    origin                 = [45.9805, -122.5030, 3350],
    ending                 = [46.3900, -122.0000, -22650],
    gridLat                = 10,
    gridLong               = 12,
    gridZ                  = 6,
    sections               = [46.12, -122.20, -1000],

    # Velocity model: 0 = 1D, 1 = 3D (file in velocityFolder)
    availableVelocity      = 1,
    namev                  = "modvMSH.txt",
    averageVelocityP       = 6,
    averageVelocityS       = 2.5,

    # Inversion
    inversionMethod        = "Tikhonov",
    MaximumIterations      = 2000,
    MaximumStallIterations = 100,
    PlotInversion          = 0,
    dampingValueQc         = [5000, 1, 1, 0.1],
    smoothingValueQc       = [0, 0, 0, 0],
    dampingValueQ          = [0.1, 0.1, 0.1, 0.1],
    smoothingValueQ        = [0, 0, 0, 0],

    # Synthetic tests: checkerboard block size and spike [lat, lon, altitude]
    sizeCheck              = 2,
    spikeLocationOrigin    = [46.05, -122.40, 0],
    spikeLocationEnd       = [46.12, -122.30, -6000],
    spikeValue             = 0.01,

    # Figures (1 = yes)
    PlotRays               = 1,
    PlotTests              = 1,
    PlotResults            = 1,
    PlotCheckers           = 1,
    PlotSpikes             = 1,
    PlotParameters         = 1,
)
