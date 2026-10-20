# Tomographic inversions (Murat_inversion, Murat_inversionQc,
# Murat_inversionQ, Murat_tikhonovQc, Murat_tikhonovQ,
# Murat_tikhonov_nonneg_fista, Murat_lsqlinQmean, Murat_smoothing,
# Murat_inputTesting, l_curve_tikh_svd).
#
# Only the "Tikhonov" inversion method is translated: damped least squares
# with non-negativity, solved by FISTA. Entry point: `murat_inversion(m)`.

"""
    tikhonov_nonneg_fista(A, b, lambda; maxit = 200, tol = 1e-6, rng) -> x

`Murat_tikhonov_nonneg_fista`: minimizes ||A x - b||² + λ²||x||² with
x ≥ 0 by FISTA. The step comes from a power-method estimate of the
Lipschitz constant started from a random vector.
This is the solver behind every Tikhonov inversion in MuRAT (Qc, Q and the
synthetic tests). MATLAB: `Murat_tikhonov_nonneg_fista.m` with the identity
regularizer (L = I), which is the only one MuRAT uses; MATLAB passes function
handles for A, A', L, L', here `A` is a matrix.
The random start uses `rng`, so results cannot be bit-identical to MATLAB
(only the step size depends on it).

# Arguments
- `A`: forward matrix (nData × nModel).
- `b`: data vector (nData).
- `lambda`: damping parameter λ (the penalty is λ²||x||²).
- `maxit`: maximum number of FISTA iterations.
- `tol`: stop when ||x_new - x|| ≤ tol · max(1, ||x||).
- `rng`: random number generator for the power-method start vector.

# Returns
- `x`: non-negative solution (nModel).
"""
function tikhonov_nonneg_fista(A, b, lambda; maxit = 200, tol = 1e-6,
                               rng = Random.MersenneTwister(0))
    n = size(A, 2)
    Ab = A' * b
    # Power method on H = A'A + λ²I: largest eigenvalue = Lipschitz constant
    # of the gradient, its inverse is the (fixed) step size
    x0 = randn(rng, n)
    x0 ./= norm(x0)
    # Hessian-vector product H v = A'A v + λ² v (matrix-free)
    H(v) = A' * (A * v) .+ lambda^2 .* v
    for _ in 1:20
        y = H(x0)
        x0 = y ./ norm(y)
    end
    Ls = max(1e-6, dot(x0, H(x0)))
    t = 1 / Ls
    # FISTA from x = 0 (feasible)
    x = zeros(n)
    yk = copy(x)
    tk = 1.0
    for _ in 1:maxit
        # Gradient step at the extrapolated point, projection onto x ≥ 0,
        # then Nesterov momentum update
        grad = H(yk) .- Ab
        xn = max.(0.0, yk .- t .* grad)
        tk1 = (1 + sqrt(1 + 4tk^2)) / 2
        yk = xn .+ ((tk - 1) / tk1) .* (xn .- x)
        if norm(xn .- x) <= tol * max(1, norm(x))
            x = xn
            break
        end
        x = xn
        tk = tk1
    end
    return x
end

"""
    l_curve_tikh_svd(U, s, d, npoints; outside_range = false) -> (rho, eta, reg_param)

L-curve of Tikhonov regularization from the SVD (Aster, Borchers and
Thurber). Returns `(rho, eta, reg_param)`: residual norms, solution norms
and regularization parameters.
Used to choose/check the damping of the Qc and Q inversions; the values are
written to `<label>/Tests/LCurve` and plotted by `murat_plot`.
MATLAB: `l_curve_tikh_svd.m` (Aster et al. toolbox, called by
`Murat_tikhonovQc.m` / `Murat_tikhonovQ.m`).

# Arguments
- `U`: left singular vectors of the forward matrix G (nData × p or nData × nData).
- `s`: singular values of G, in decreasing order (length p).
- `d`: data vector (nData).
- `npoints`: number of regularization parameters on the curve.
- `outside_range`: add the part of `d` outside the range of `U` to the
  residual norm (off by default, as in MuRAT, see the comment below).

# Returns
- `rho`: residual norms ||G m - d|| (npoints).
- `eta`: solution norms ||m|| (npoints).
- `reg_param`: regularization parameters, from s[1] down to
  max(s[p], 16 eps s[1]), logarithmically spaced (npoints).
"""
function l_curve_tikh_svd(U, s, d, npoints; outside_range = false)
    m, n = size(U)
    p = length(s)
    # Data projected on the left singular vectors
    dp = U' * d
    # MuRAT passes the full U of svd(G), so this term is ~0 there: it is
    # only added if outside_range = true
    dr = outside_range ? norm(d)^2 - norm(dp)^2 : 0.0
    dp = dp[1:p]
    dps = dp ./ s
    # Regularization parameters: geometric sequence from s[1] (i = 1) to the
    # smallest singular value (i = npoints)
    reg = zeros(npoints)
    reg[npoints] = max(s[p], s[1] * 16 * eps())
    ratio = (s[1] / reg[npoints])^(1 / (npoints - 1))
    for i in npoints-1:-1:1
        reg[i] = ratio * reg[i+1]
    end
    eta = zeros(npoints)
    rho = zeros(npoints)
    s2 = s .^ 2
    # Tikhonov filter factors f = s²/(s² + α²): solution norm from the
    # filtered coefficients, residual norm from the filtered-out part
    for i in 1:npoints
        f = s2 ./ (s2 .+ reg[i]^2)
        eta[i] = norm(f .* dps)
        rho[i] = norm((1 .- f) .* dp)
    end
    (m > n && dr > 0) && (rho = sqrt.(rho .^ 2 .+ dr))
    return rho, eta, reg
end

"""
    tikhonov_solve(G, b, dampValue, x0, maxIt; rng) -> (x, misfit, lcurve)

Common Tikhonov step of the Qc and Q inversions: L-curve data, non-negative
FISTA solution and relative misfit.
MATLAB: `Murat_tikhonovQc.m` / `Murat_tikhonovQ.m` (without the L-curve
figure, which `murat_plot` draws from `lcurve`; the economy SVD gives the same
L-curve as MATLAB's full SVD because the residual outside the range is not
added). The misfit is computed as in MATLAB: norm of the squared residuals.

# Arguments
- `G`: forward matrix (nData × nModel), already weighted for Qc.
- `b`: data vector (nData), weighted like `G`.
- `dampValue`: damping parameter λ.
- `x0`: starting model (nModel), used only as reference for the misfit.
- `maxIt`: maximum FISTA iterations.
- `rng`: random number generator passed to the FISTA solver.

# Returns
- `x`: non-negative solution (nModel).
- `misfit`: relative misfit ||(b - G x).²|| / ||(b - G x0).²|| (< 1 means the
  solution fits better than the starting model).
- `lcurve`: `(rho, eta, reg_param)` from `l_curve_tikh_svd` (100 points).
"""
function tikhonov_solve(G, b, dampValue, x0, maxIt; rng)
    F = svd(G, full = false)
    lcurve = l_curve_tikh_svd(F.U, F.S, b, 100)
    obj0 = norm((b .- G * x0) .^ 2)
    x = tikhonov_nonneg_fista(G, b, dampValue; maxit = maxIt, tol = 1e-7, rng = rng)
    obj = norm((b .- G * x) .^ 2)
    return x, obj / obj0, lcurve
end

"""
    murat_smoothing(coordPrior, smoothValue, sigmaM) -> C

`Murat_smoothing`: Gaussian smoothing (covariance) matrix between grid
nodes (`coordPrior`: lat, lon in degrees, z in m).
MATLAB: `Murat_smoothing.m`. The distance between nodes is the great-circle
(haversine) horizontal distance combined with the vertical difference.
Deliberate difference: when `smoothValue === nothing`, MATLAB indexes
`nearestVals(idx(1))`, which can fail or pick the wrong node; here the first
node at the nearest different depth is used.
Note: in the Tikhonov method the matrix is computed by `murat_inversion_qc`
but not used by the solver (as in MATLAB, where the smoothing term is
commented out in `Murat_tikhonovQc.m`).

# Arguments
- `coordPrior`: N × 3 node coordinates; columns 1-2 in degrees (read as lat,
  lon as documented in MATLAB), column 3 in m.
- `smoothValue`: correlation length in km, or `nothing` for 3 × the depth step
  of the grid.
- `sigmaM`: model variance (scales the whole matrix).

# Returns
- `C`: N × N matrix sigmaM · exp(-D²/(4 smoothValue²)), D = node distance in
  km; sparse (entries with D > 3 smoothValue set to zero) when N > 2000.
"""
function murat_smoothing(coordPrior, smoothValue, sigmaM)
    # Mean Earth radius (km)
    Rk = 6371.0
    lat = deg2rad.(coordPrior[:, 1])
    lon = deg2rad.(coordPrior[:, 2])
    z = coordPrior[:, 3] ./ 1000
    # Pairwise horizontal distances (km) by the haversine formula
    dlat = lat .- lat'
    dlon = lon .- lon'
    aa = sin.(dlat ./ 2) .^ 2 .+ (cos.(lat) .* cos.(lat')) .* sin.(dlon ./ 2) .^ 2
    dh = Rk .* 2 .* atan.(sqrt.(aa), sqrt.(1 .- aa))
    # 3D straight-line distances (km)
    D = sqrt.(dh .^ 2 .+ (z .- z') .^ 2)
    # Default correlation length: 3 × distance (km) from the first node's depth
    # to the nearest different depth
    if smoothValue === nothing
        v = coordPrior[:, 3]
        d = abs.(v .- v[1])
        d[v .== v[1]] .= Inf
        smoothValue = abs(3 * (v[findfirst(==(minimum(d)), d)] - v[1])) / 1000
    end
    smoothValue = abs(smoothValue)
    # Gaussian covariance; exp(-0.5 D²/(2 s²)) is kept as written in MATLAB
    C = sigmaM .* exp.(-0.5 .* D .^ 2 ./ (2 * smoothValue^2))
    size(coordPrior, 1) > 2000 && (C = sparse((D .<= 3smoothValue) .* C))
    return C
end

"""
    murat_inversion_qc(Ac, Qm, method, iter, RZZ, coordPrior, dampValue, smoothValue; rng)
        -> (Qc, misfit, dampValue, lcurve)

`Murat_inversionQc` (Tikhonov): weighted, damped, non-negative inversion of
the measured inverse Qc. Returns `(Qc, misfit, dampValue, lcurve)`.
Maps the inverse coda quality factor measured on each trace onto the grid
cells, using the coda sensitivity kernels (Qc⁻¹ of a trace ≈ kernel-weighted
average of the cell values). Data are weighted by 1/uncertainty.
MATLAB: `Murat_inversionQc.m`; only `method == "Tikhonov"` is available
(`"Particle"`, `"Annealing"`, `"Global"` raise an error). MATLAB's lower bound
min(Qm) applies only to those other methods; FISTA uses x ≥ 0.

# Arguments
- `Ac`: Qc inversion matrix, rows = retained traces, columns = crossed cells,
  each row a normalized kernel summing to 1 (nData × nModel).
- `Qm`: measured inverse Qc per trace (nData).
- `method`: inversion method name (`s.inversionMethod`).
- `iter`: maximum solver iterations.
- `RZZ`: uncertainty of each Qm measurement (nData); weights are 1/RZZ.
- `coordPrior`: coordinates of the inverted cells (nModel × 3), only used for
  the smoothing matrix.
- `dampValue`: damping λ, or `nothing` for mean(m0)/10.
- `smoothValue`: smoothing length (km), `nothing` or 0 for no smoothing matrix.
- `rng`: random number generator for the FISTA solver.

# Returns
- `Qc`: inverse Qc per cell (nModel).
- `misfit`: relative misfit reduction with respect to the starting model.
- `dampValue`: damping actually used (needed by the synthetic tests).
- `lcurve`: `(rho, eta, reg_param)`.
"""
function murat_inversion_qc(Ac, Qm, method, iter, RZZ, coordPrior, dampValue, smoothValue; rng)
    # Data weighting: W = diag(1/RZZ²) multiplies both G and the data
    w = 1 ./ RZZ
    W = Diagonal(w .^ 2)
    # Starting model (regionalization): kernel-weighted mean of the data per cell
    m0 = vec(sum(Ac .* Qm, dims = 1)) ./ vec(sum(Ac, dims = 1))
    # Smoothing covariance with model variance (mean(m0)/10)²; computed as in
    # MATLAB but not used by the Tikhonov solver
    CovM = (smoothValue !== nothing && smoothValue != 0) ?
           murat_smoothing(coordPrior, smoothValue, mean(m0 ./ 10)^2) : 0
    # Default damping: 1/10 of the mean starting inverse Qc
    dampValue = dampValue === nothing ? mean(m0 ./ 10) : abs(dampValue)
    method == "Tikhonov" ||
        error("Inversion method '$method' is not available yet in the Julia version (use \"Tikhonov\").")
    x, misfit, lcurve = tikhonov_solve(W * Ac, W * Qm, dampValue, m0, iter; rng = rng)
    return x, misfit, dampValue, lcurve
end

"""
    murat_lsqlin_qmean(cf, l, time0, rapsp, Qc, tCm, te) -> (d0, Q0, mDest)

`Murat_lsqlinQmean`: average coda constant and Q from the coda-normalized
energy ratios, data for the Q inversion and estimated diffusion constant.
Returns `(d0, Q0, mDest)`, with `Q0 = [estimates variances]` (2×2).
Coda normalization: log(E_direct/E_coda · l²)/(2π f) = c - t Q⁻¹, fitted by
least squares over all traces (c = coda constant, Q⁻¹ = average inverse Q).
The diffusion constant D is then searched on a grid so that the diffusive
coda-energy model reproduces the estimated c (best D per trace, averaged).
MATLAB: `Murat_lsqlinQmean.m` (unconstrained `lsqlin` = least squares here).

# Arguments
- `cf`: central frequency (Hz).
- `l`: total ray length per trace (km).
- `time0`: travel time per trace (s).
- `rapsp`: direct-to-coda energy ratio per trace.
- `Qc`: inverse Qc per trace.
- `tCm`: coda start lapse time per trace (s).
- `te`: coda end lapse time per trace (s), `tCm + codaWindow`.

# Returns
- `d0`: data for the Q inversion, log term minus the coda constant (nData).
- `Q0`: 2 × 2; row 1 = coda constant, row 2 = average inverse Q;
  column 1 = estimate, column 2 = variance.
- `mDest`: mean estimated diffusion constant (km²/s).
"""
function murat_lsqlin_qmean(cf, l, time0, rapsp, Qc, tCm, te)
    # Linear model data = c - t Q⁻¹ with unknowns [c, Q⁻¹]
    data = log.(rapsp .* l .^ 2) ./ (2π * cf)
    G = hcat(ones(length(time0)), -time0)
    x = G \ data
    # Parameter variances from the residual variance (dof = nData - 2)
    res = G * x .- data
    dof = max(size(G, 1) - 2, 1)
    cov_x = (res' * res) / dof .* pinv(G' * G)
    Q0 = hcat(x, diag(cov_x))
    d0 = data .- Q0[1, 1]
    cEst = Q0[1, 1]
    # Grid search of D (linear grid 1e-3 to 1e3, 2000 values): for each trace
    # the D whose diffusive coda constant (1.5 log(4πD) - log A)/(2πf) is
    # closest to cEst; A = difference of the diffusion energy at the two
    # coda times. MATLAB writes (tCm + te)^(-1.5) in the first term: kept.
    # Ties keep the first (smallest) D and traces with A ≤ 0 everywhere get
    # Ds[1], as MATLAB's min.
    Ds = collect(range(1e-3, 1e3, length = 2000))
    Dest = zeros(length(l))
    for n in eachindex(l)
        best, bestr = Ds[1], Inf
        for D in Ds
            A = (tCm[n] + te[n])^(-1.5) * exp(-l[n]^2 / (4D * te[n]) - 2π * cf * Qc[n] * te[n]) -
                tCm[n]^(-1.5) * exp(-l[n]^2 / (4D * tCm[n]) - 2π * cf * Qc[n] * tCm[n])
            A > 0 || continue
            r = abs(cEst - (1.5 / (2π * cf) * log(4π * D) - 1 / (2π * cf) * log(A)))
            if isfinite(r) && r < bestr
                best, bestr = D, r
            end
        end
        Dest[n] = best
    end
    return d0, Q0, mean(Dest)
end

"""
    murat_inversion_q(A, d0, rapsp, method, Q0, iter, coordPrior, dampValue, smoothValue; rng)
        -> (Q, misfit, dampValue, lcurve)

`Murat_inversionQ` (Tikhonov): inversion of the coda-normalized energy
ratios for Q along the rays. Returns `(Q, misfit, dampValue, lcurve)`.
The data `d0` (from `murat_lsqlin_qmean`, ≈ -t Q⁻¹ per trace) are modelled
as `A * Qinv`, with A = minus the travel time spent in each cell, so the
result is the total-attenuation (inverse Q) value per cell.
MATLAB: `Murat_inversionQ.m`; only `method == "Tikhonov"` is available. In
the Tikhonov branch MATLAB does not use `Q0`, `coordPrior` and `smoothValue`
(they serve the other methods / smoothing), so they are unused here too.

# Arguments
- `A`: Q matrix, minus the travel time (s) of each ray in each cell; rows =
  traces, columns = crossed cells (nData × nModel).
- `d0`: coda-normalized data (nData).
- `rapsp`: direct-to-coda energy ratios (nData), used for the default damping.
- `method`: inversion method name (`s.inversionMethod`).
- `Q0`: average coda constant and inverse Q (2 × 2, see `murat_lsqlin_qmean`).
- `iter`: maximum solver iterations.
- `coordPrior`: coordinates of the inverted cells (unused by Tikhonov).
- `dampValue`: damping λ, or `nothing` for 0.1 × mean(rapsp).
- `smoothValue`: smoothing length (unused by Tikhonov).
- `rng`: random number generator for the FISTA solver.

# Returns
- `Q`: inverse Q per cell (nModel).
- `misfit`: relative misfit reduction with respect to the starting model.
- `dampValue`: damping actually used.
- `lcurve`: `(rho, eta, reg_param)`.
"""
function murat_inversion_q(A, d0, rapsp, method, Q0, iter, coordPrior, dampValue, smoothValue; rng)
    MQ = size(A, 2)
    # Starting model (regionalization): d = d0 / (-travel time) = Q⁻¹ per
    # trace, averaged per cell with weights = fraction of the cell's total time
    d = d0 ./ vec(sum(A, dims = 2))
    Ab = abs.(A ./ sum(A, dims = 1))
    m0 = [sum(Ab[Ab[:, i] .> 0, i] .* d[Ab[:, i] .> 0]) for i in 1:MQ]
    # Default damping: 0.1 × mean energy ratio (as MATLAB)
    dampValue = dampValue === nothing ? 0.1 * mean(rapsp) : abs(dampValue)
    method == "Tikhonov" ||
        error("Inversion method '$method' is not available yet in the Julia version (use \"Tikhonov\").")
    x, misfit, lcurve = tikhonov_solve(A, d0, dampValue, m0, iter; rng = rng)
    return x, misfit, dampValue, lcurve
end

"""
    murat_input_testing(I, spike_o, spike_e, x, y, z) -> (checkInput, spike)

`Murat_inputTesting`: checkerboard pattern in MuRAT list order and cells of
the spike (between `spike_o` and `spike_e`, given as lat, lon, altitude).
MATLAB: `Murat_inputTesting.m`.

# Arguments
- `I`: checkerboard array (nx × ny × nz) of 0/1 from `checkerboard3d`.
- `spike_o`, `spike_e`: opposite corners of the spike box, `[lat, lon,
  altitude (m)]` (`spike_o` is the top), or `spike_o === nothing` for no spike.
- `x`, `y`, `z`: grid vectors (longitude °, latitude °, altitude m).

# Returns
- `checkInput`: checkerboard values per cell in MuRAT list order (x slowest,
  z fastest).
- `spike`: Bool mask of the cells strictly inside the box, or `nothing`.
"""
function murat_input_testing(I, spike_o, spike_e, x, y, z)
    # (nx, ny, nz) array to list with z fastest, as reshape(permute(I,[3 2 1]),[],1)
    checkInput = vec(permutedims(I, (3, 2, 1)))
    spike_o === nothing && return checkInput, nothing
    # r columns: lon, lat, altitude, compared with spike [lat, lon, altitude]
    r = unfold_xyz(x, y, z)
    spike = (r[:, 1] .> spike_o[2]) .& (r[:, 1] .< spike_e[2]) .&
            (r[:, 2] .> spike_o[1]) .& (r[:, 2] .< spike_e[1]) .&
            (r[:, 3] .< spike_o[3]) .& (r[:, 3] .> spike_e[3])
    return checkInput, spike
end

"""
    valueat(v, k)

`v[k]`, or `nothing` if `v === nothing`: per-frequency optional input values
(damping, smoothing), as MATLAB's `if ~isempty(v), v(k), else [] end`.
No MATLAB counterpart (Julia helper).
"""
valueat(v, k) = v === nothing ? nothing : v[k]

"""
    murat_inversion(m::Murat; seed = 0) -> Murat

`Murat_inversion`: peak-delay regionalization, Qc and Q inversions,
checkerboard and spike tests, for each frequency. Writes the result tables
in `<label>/TXT` and the L-curve data in `<label>/Tests/LCurve`.
Third stage of a run, after `murat_data` (measurements, rays, kernels,
selection) and before `murat_summary` / `murat_plot`.
MATLAB: `Murat_inversion.m`. Differences: the L-curves are stored in
`m.Qc.lcurves` and written as text (MATLAB saves figures); the solver
outputs/exit flags of the other methods are not stored; the FISTA random
start comes from `seed`.

# Arguments
- `m`: run after `murat_data`, modified in place.
- `seed`: seed of the random generator used by the FISTA solver.

# Returns
- `m`, with `m.PD.modvPeakDelay` (nCells × 5 × nFreq), `m.Qc.modvQc` and
  `m.Q.modvQ` (nCells × 9 × nFreq), relative misfits, damping values,
  L-curves, coda constants, average Q and diffusion constant. Columns of
  `modv*[:, :, f]`: 1-3 lon, lat, altitude (m); 4 model (peak-delay
  variation, inverse Qc or inverse Q); 5 hit count (peak delay only, zero
  for Qc and Q); 6-7 checkerboard input/output; 8-9 spike input/output.
"""
function murat_inversion(m::Murat; seed = 0)
    s = m.input
    g = m.grid
    rng = Random.MersenneTwister(seed)
    cf = s.centralFrequency
    modvP = g.modvPlot
    x, y, z = g.x, g.y, g.z
    nxyz = length(x) * length(y) * length(z)
    tWm = s.codaWindow
    # Diffusion constant D = v_S / (3 Le⁻¹ B0) per frequency (Wu 1985; Wu and
    # Aki 1988), from the S velocity, inverse extinction length and albedo
    D = s.averageVelocityS ./ 3 ./ s.iExtinctionLength ./ s.albedo

    # Data, matrices and selection masks from `murat_data` (rows = traces,
    # columns of the matrices = grid cells; rcXX = cells crossed, retainXX =
    # traces kept, both per frequency)
    Apd_i = m.PD.inversionMatrixPeakDelay
    Ac_i = m.Qc.inversionMatrixQc
    A_i = m.Q.inversionMatrixQ
    luntot = m.rays.totalLengthRay
    time0 = m.rays.travelTime
    Qm = m.Qc.inverseQc
    RZZ = m.Qc.uncertaintyQc
    lpdelta = m.PD.variationPeakDelay
    rapE = m.Q.energyRatioBodyCoda
    retainPD, retainQc, retainQ = m.PD.retainPeakDelay, m.Qc.retainQc, m.Q.retainQ
    rcPD, rcQc, rcQ = m.PD.raysPeakDelay, m.Qc.raysQc, m.Q.raysQ
    tCoda = m.Qc.tCoda
    nFreq = size(rcPD, 2)

    # Output models: modvPD columns [lon lat alt ΔlogPD hitCount]; modvQc and
    # modvQ columns [lon lat alt model hitC(0) checkIn checkOut spikeIn spikeOut]
    modvPD = zeros(nxyz, 5, nFreq)
    modvQc = zeros(nxyz, 9, nFreq)
    modvQ = zeros(nxyz, 9, nFreq)
    residualQc = zeros(nFreq)
    residualQ = zeros(nFreq)
    dValueQc = zeros(nFreq)
    dValueQ = zeros(nFreq)
    cConstants = zeros(2, nFreq)
    lcurves = Dict{String,Any}()
    averageQ = nothing
    averageD = nothing

    # Checkerboard (0/1 blocks of sizeCheck cells) and spike mask, frequency
    # independent (MATLAB builds them inside the loop, once)
    I = checkerboard3d((length(x), length(y), length(z)), s.sizeCheck)
    checkInput, spikeInput = murat_input_testing(I, s.spikeLocationOrigin, s.spikeLocationEnd, x, y, z)

    for k in 1:nFreq
        cf_k = cf[k]
        fcName = replace(string(fmt_g(cf_k)), "." => "_")
        D_k = D[k]
        modvPD[:, 1:3, k] .= modvP[:, 1:3]
        modvQc[:, 1:3, k] .= modvP[:, 1:3]
        modvQ[:, 1:3, k] .= modvP[:, 1:3]

        # Peak-delay regionalization
        # (no inversion: each cell gets the mean log10 peak-delay variation of
        # the rays crossing it; cells without data get the mean of all cells)
        rcpd, rtpd = rcPD[:, k], retainPD[:, k]
        Apd_k = Apd_i[rtpd, rcpd]
        boxes = Apd_k .> 0
        counts = vec(sum(boxes, dims = 1))
        mpd = vec(sum(boxes .* lpdelta[rtpd, k], dims = 1)) ./ counts
        mpd[isnan.(mpd)] .= mean(filter(!isnan, mpd))
        modvPD[rcpd, 4, k] .= mpd
        modvPD[rcpd, 5, k] .= counts

        # Qc inversion
        rcqc, rtqc = rcQc[:, k], retainQc[:, k]
        Ac_k = Ac_i[rtqc, rcqc, k]
        RZZ_k = RZZ[rtqc, k]
        Qm_k = Qm[rtqc, k]
        coordQc = modvQc[rcqc, 1:3, k]
        sol, fval, dQc, lc = murat_inversion_qc(Ac_k, Qm_k, s.inversionMethod, s.MaximumIterations,
                                                RZZ_k, coordQc, valueat(s.dampingValueQc, k),
                                                valueat(s.smoothingValueQc, k); rng = rng)
        modvQc[rcqc, 4, k] .= sol
        residualQc[k] = fval
        lcurves["Qc_$(fcName)_Hz"] = lc
        @printf("Qc relative misfit reduction at %gHz is: %.4f\n", cf_k, fval)

        # Q inversion
        # (traces retained for both Q and Qc, because coda normalization uses Qc)
        rcq, rtq = rcQ[:, k], retainQ[:, k]
        keep = rtq .& rtqc
        A_k = A_i[keep, rcq]
        Qc_k = Qm[keep, k]
        l = luntot[keep] ./ 1000
        time0_k = time0[keep]
        rapE_k = rapE[keep, k]
        tCm = tCoda[keep, k]
        coordQ = modvQ[rcq, 1:3, k]
        te = tCm .+ tWm
        # A priori coda constant per trace from the diffusion model with D_k
        # (only printed and stored for comparison); Complex because the log
        # argument can be negative, as MATLAB's log; the real part is used
        cCoda = 1.5 / 2 / π / cf_k * log(4π * D_k) .-
                1 / 2 / π / cf_k .* log.(Complex.(te .^ (-1.5) .* exp.(-l .^ 2 ./ 4 ./ D_k ./ te .-
                                                                   2π * cf_k .* Qc_k .* te) .-
                                                tCm .^ (-1.5) .* exp.(-l .^ 2 ./ 4 ./ D_k ./ tCm .-
                                                                      2π * cf_k .* Qc_k .* tCm)))
        # Estimated coda constant, average Q⁻¹ and diffusion constant
        d0, Q0, mDest = murat_lsqlin_qmean(cf_k, l, time0_k, rapE_k, Qc_k, tCm, te)
        cEst = Q0[1, 1]
        cConstants[:, k] .= (mean(real.(cCoda)), cEst)
        @printf("A priori and estimated coda constants are %.4g and %.4g at %g Hz\n",
                mean(real.(cCoda)), cEst, cf_k)
        @printf("A priori and estimated diffusion constants are %.4g and %.4g at %g Hz\n",
                D_k, mDest, cf_k)
        # Only the last frequency's values are kept (as MATLAB)
        averageQ, averageD = Q0, mDest
        solQ, fvalQ, dQ, lcQ = murat_inversion_q(A_k, d0, rapE_k, s.inversionMethod, Q0,
                                                 s.MaximumIterations, coordQ,
                                                 valueat(s.dampingValueQ, k),
                                                 valueat(s.smoothingValueQ, k); rng = rng)
        modvQ[rcq, 4, k] .= solQ
        residualQ[k] = fvalQ
        lcurves["Q_$(fcName)_Hz"] = lcQ
        @printf("Q relative misfit reduction at %gHz is: %.4f\n", cf_k, fvalQ)

        # Checkerboard
        # Input: checkerboard alternating min and max measured Qc⁻¹ (col 6);
        # spike input = mean Qc⁻¹ background with spikeValue in the box (col 8)
        modvQc[checkInput .== 1, 6, k] .= minimum(Qm_k)
        modvQc[checkInput .== 0, 6, k] .= maximum(Qm_k)
        modvQc[:, 8, k] .= mean(Qm_k)
        spikeInput === nothing || (modvQc[spikeInput, 8, k] .= s.spikeValue)

        # Synthetic Qc data from the input, inverted with the same matrix,
        # weights and damping as the real data (col 7 = recovered model)
        re_check = Ac_k * modvQc[rcqc, 6, k]
        solc, _, _, _ = murat_inversion_qc(Ac_k, re_check, s.inversionMethod, s.MaximumIterations,
                                           RZZ_k, coordQc, dQc, valueat(s.smoothingValueQc, k); rng = rng)
        modvQc[rcqc, 7, k] .= solc

        # Q tests reuse the Qc inputs (cols 6-8). Synthetic energy ratios
        # E = exp(2πf (c + ∫Q⁻¹)) / l² go through the same coda-normalization
        # fit and Q inversion as the real data
        modvQ[:, 6:8, k] .= modvQc[:, 6:8, k]
        re_checkQ = A_k * modvQ[rcq, 6, k]
        synth = exp.(2π * cf_k .* (cEst .+ re_checkQ)) ./ l .^ 2
        d0c, Q0c, _ = murat_lsqlin_qmean(cf_k, l, time0_k, synth, Qc_k, tCm, te)
        solcQ, _, _, _ = murat_inversion_q(A_k, d0c, rapE_k, s.inversionMethod, Q0c,
                                           s.MaximumIterations, coordQ, dQ,
                                           valueat(s.smoothingValueQ, k); rng = rng)
        modvQ[rcq, 7, k] .= solcQ

        # Spike
        # (same procedure as the checkerboard, input col 8, output col 9)
        if spikeInput !== nothing
            re_spike = Ac_k * modvQc[rcqc, 8, k]
            sols, _, _, _ = murat_inversion_qc(Ac_k, re_spike, s.inversionMethod, s.MaximumIterations,
                                               RZZ_k, coordQc, dQc, valueat(s.smoothingValueQc, k); rng = rng)
            modvQc[rcqc, 9, k] .= sols
            re_spikeQ = A_k * modvQ[rcq, 8, k]
            synth = exp.(2π * cf_k .* (cEst .+ re_spikeQ)) ./ l .^ 2
            d0s, Q0s, _ = murat_lsqlin_qmean(cf_k, l, time0_k, synth, Qc_k, tCm, te)
            solsQ, _, _, _ = murat_inversion_q(A_k, d0s, rapE_k, s.inversionMethod, Q0s,
                                               s.MaximumIterations, coordQ, dQ,
                                               valueat(s.smoothingValueQ, k); rng = rng)
            modvQ[rcq, 9, k] .= solsQ
        end

        dValueQc[k] = dQc
        dValueQ[k] = dQ

        # Result tables (same names and columns as MATLAB) and L-curve values
        write_table(outdir(m, "TXT", "peakdelay_$(fcName)_Hz.txt"),
                    ["Lon", "Lat", "Depth", "ΔlogPD", "HitC"], modvPD[:, :, k])
        hdr9 = ["Check In", "Check Out", "Spike In", "Spike Out"]
        write_table(outdir(m, "TXT", "Qc_$(fcName)_Hz.txt"),
                    vcat(["Lon", "Lat", "Depth", "Qc", "HitC"], hdr9), modvQc[:, :, k])
        write_table(outdir(m, "TXT", "Q_$(fcName)_Hz.txt"),
                    vcat(["Lon", "Lat", "Depth", "Q", "HitC"], hdr9), modvQ[:, :, k])
        for (name, lc) in (("Qc", lc), ("Q", lcQ))
            write_table(outdir(m, "Tests", "LCurve", "L-curve_$(name)_$(fcName)_Hz.txt"),
                        ["rho", "eta", "reg_param"], hcat(lc...))
        end
    end

    # Store the results in the run (same fields as the MATLAB structure;
    # DiffConstant goes to m.grid instead of Murat.input)
    g.DiffConstant = D
    m.PD.modvPeakDelay = modvPD
    m.Qc.modvQc = modvQc
    m.Qc.residualRedQc = residualQc
    m.Qc.dampingValueQc = dValueQc
    m.Qc.lcurves = lcurves
    m.Q.residualRedQ = residualQ
    m.Q.dampingValueQ = dValueQ
    m.Q.modvQ = modvQ
    m.Q.codaConstants = cConstants
    m.Q.averageQ = averageQ
    m.Q.estimatedDiffusion = averageD
    return m
end

"""
    fmt_g(x::Real)

`x` as an `Int` if it is integer-valued, so that frequencies print like
MATLAB `num2str` in file names (`6` rather than `6.0`, `1.5` stays `1.5`).
No MATLAB counterpart (Julia helper).
"""
fmt_g(x::Real) = isinteger(x) ? Int(x) : x

"""
    write_table(path, header, M)

Tab-delimited table with a header line (as MATLAB `writecell`).
Writes `header` (vector of column names) and the rows of the numeric
matrix `M`, each value with 15 significant digits. Used for the `TXT` result
tables and the L-curve files. MATLAB: `writecell(..., 'Delimiter', '\\t')` in
`Murat_inversion.m`.
"""
function write_table(path, header, M)
    open(path, "w") do io
        println(io, join(header, '\t'))
        for r in axes(M, 1)
            println(io, join([@sprintf("%.15g", v) for v in M[r, :]], '\t'))
        end
    end
end
