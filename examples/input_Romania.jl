# Romania example (translated from Murat_inputRomania.m of MuRAT 4.0); needs `using NCDatasets` (LITHO1.0)
# Same field names as Murat.input.<name> in MATLAB
# (`end` is called `ending`; MATLAB [] is `nothing`).
# The data are in this folder (examples/); the results go to the folder from
# which the run is started.

using MuRAT

MuratInput(
    # Data folder and output label
    dataDirectory          = joinpath(@__DIR__, "sac_Romania"),
    label                  = "Romania",
    velocityFolder         = joinpath(@__DIR__, "velocity_models"),

    # SAC header times: origin (nothing = "o"), P and S picks (nothing = "t0")
    originTime             = "o",
    PTime                  = "t0",
    STime                  = "t0",

    # Phase used for picking: 2 = P, 3 = S
    POrS                   = 3,

    # Central frequencies (Hz) and number of components
    centralFrequency       = [3, 6, 12, 18],
    envelopeSmoothTime     = 1,
    components             = 1,

    # Declustering factor (nothing = no declustering)
    declustering           = 5,

    # Peak delay: limits (s) and outlier threshold (standard deviations)
    minimumPeakDelay       = 0.1,
    maximumPeakDelay       = 10,
    stDevPD                = 2,

    # Coda normalization: direct-wave window, noise start (s), noise threshold
    bodyWindow             = 1,
    startNoise             = 5,
    tresholdNoise          = 3,

    # Coda: lapse-time method ("Constant", "Peak" or "Travel"), start, window
    lapseTimeMethod        = "Constant",
    startLapseTime         = 30,
    codaWindow             = 30,
    mintravel              = 0.2,
    maxtravel              = Inf,

    # Kernels: albedo and inverse extinction length per frequency
    albedo                 = [0.5, 0.5, 0.5, 0.5],
    iExtinctionLength      = [0.02, 0.02, 0.02, 0.02],
    kernelTreshold         = 2,
    # "differential" (MuRAT 4.0, default) or "middle" (MuRAT 3)
    kernelType             = "differential",

    # Qc measurement ("Linearized" or "NonLinear") and fit threshold
    QcMeasurement          = "Linearized",
    fitTresholdLinear      = 0.3,

    # Grid: origin and end [lat, lon, altitude (m)], number of nodes
    origin                 = [43, 20, 1000],
    ending                 = [49, 29, -50000],
    gridLat                = 5,
    gridLong               = 9,
    gridZ                  = 6,
    sections               = [45.5, 26.7, -7000],

    # Velocity model: 0 = 1D, 1 = 3D (file in velocityFolder)
    availableVelocity      = 0,
    namev                  = "LITHO1.0.nc",
    averageVelocityP       = 6,
    averageVelocityS       = 3,

    # Inversion
    inversionMethod        = "Tikhonov",
    MaximumIterations      = 2000,
    MaximumStallIterations = 100,
    PlotInversion          = 0,
    dampingValueQc         = [0.1, 0.1, 0.1, 0.1],
    smoothingValueQc       = nothing,
    dampingValueQ          = [0.01, 0.01, 0.01, 0.01],
    smoothingValueQ        = nothing,

    # Synthetic tests: checkerboard block size and spike [lat, lon, altitude]
    sizeCheck              = 2,
    spikeLocationOrigin    = [44, 23, -10000],
    spikeLocationEnd       = [46, 25, -30000],
    spikeValue             = 0.01,

    # Figures (1 = yes)
    PlotRays               = 1,
    PlotTests              = 1,
    PlotResults            = 1,
    PlotCheckers           = 1,
    PlotSpikes             = 1,
    PlotParameters         = 1,
)
