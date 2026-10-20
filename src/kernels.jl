# Coda sensitivity kernels (Murat_paasschensFunction, Murat_kernels,
# Murat_codaMatrix). Original creators: Angel De La Torre & Edoardo Del Pezzo.
#
# Units in this file: distances r, R, DR, D0 in km; coordinates of grids, events
# and stations in m; velocity v in km/s; times in s; inverse extinction length
# Le_1 in 1/km; albedo B0 dimensionless.

"""
    paasschens(r, v, B0, Le_1, dt, T) -> (t0, A0, N, coda)

Paasschens approximation of the energy Green's function (radiative transfer with
isotropic scattering) at distance `r` from the source: a delta (direct, coherent
energy) of amplitude `A0` at `t0 = r/v`, plus a coda sampled every `dt` from `t0`
to `T + 2dt`. Building block of the coda sensitivity kernel (`murat_kernels`).
MATLAB: `Murat_paasschensFunction.m` (the time vector output `t` is not returned).

# Arguments
- `r`: source distance (km).
- `v`: average velocity (km/s).
- `B0`, `Le_1`: albedo and inverse extinction length (1/km), vectors (one value per
  frequency).
- `dt`: time step (s).
- `T`: lapse time of interest (s).

# Returns
- `t0`: arrival time of the delta = start of the coda (s).
- `A0`: delta amplitude, `exp(-Le_1 v t0) / (4π r² v)`, one value per frequency.
- `N`: number of coda samples.
- `coda`: N × length(B0) coda samples at `t0, t0 + dt, ...`.
"""
function paasschens(r, v, B0, Le_1, dt, T)
    t0 = r / v
    t = matlab_colon(t0, dt, T + 2dt)
    N = length(t)
    A0 = exp.(-Le_1 .* v .* t0) ./ (4π * r^2 * v)
    lR = length(B0)
    coda = zeros(N, lR)
    # Coda E(r, t) = f1 * f2 * f3 * F, with F(x) = exp(x) sqrt(1 + 2.026/x)
    # (Paasschens' approximation); one_minus = 1 - t0²/t² clipped at 0
    for c in 1:lR, n in 1:N
        one_minus = max(1 - t0^2 / t[n]^2, 0.0)
        f1 = one_minus^(1 / 8)
        f2 = (3 * B0[c] * Le_1[c] / (4π * v) / t[n])^(3 / 2)
        f3 = exp(-Le_1[c] * v * t[n])
        x = v * t[n] * B0[c] * Le_1[c] * one_minus^(3 / 4)
        F = x > 1e-30 ? exp(x) * sqrt(1 + 2.026 / x) : 0.0
        coda[n, c] = f1 * f2 * f3 * F
    end
    # The sample at t0 is singular but integrable: MuRAT stores 2/3 of the next
    # sample there (MATLAB `coda(1)` is linear indexing, hence first column only)
    N >= 2 && (coda[1, 1] = coda[2, 1] * 2 / 3)   # as in MuRAT: first column only
    return t0, A0, N, coda
end

"""
    murat_kernels(T, event, station, modv, v, kT, B0, Le_1, lapseTimeMethod) -> (K, r)

Coda sensitivity kernel at lapse time `T` for an event-station couple, on a
regular grid around them (spacing = inversion grid step / kT). The kernel at a
point at distances r1 (from the source) and r2 (from the station) is the
convolution of the two Paasschens energy Green's functions evaluated at `T`; it is
non-zero only inside the ellipsoid r1 + r2 < v T. Used in `murat_data` to build the
Qc inversion matrix (through `murat_coda_matrix`).
MATLAB: `Murat_kernels.m`. Differences: no interactive pause when the grid has more
than 50e6 points (only a warning); K is points × frequencies (MATLAB points × 1 ×
frequencies).

# Arguments
- `T`: lapse time (s); a vector (one per frequency) is averaged when
  `lapseTimeMethod == "Peak"`.
- `event`, `station`: coordinates `[x, y, z]` in m (only the first 3 are used).
- `modv`: inversion grid, nModel × ≥3 `[x y z ...]` in m (sets the grid step).
- `v`: average (S-wave) velocity (km/s).
- `kT`: factor dividing the inversion grid step (`kernelTreshold` setting); larger
  is finer and slower.
- `B0`, `Le_1`: albedo and inverse extinction length (1/km), one per frequency.
- `lapseTimeMethod`: `"Peak"`, `"Constant"` or `"Travel"`.

# Returns
- `K`: Nxyz × nFreq kernel values (zero outside the ellipsoid).
- `r`: Nxyz × 3 kernel grid coordinates in m (same reference as `event`), MuRAT list
  order.
"""
function murat_kernels(T, event, station, modv, v, kT, B0, Le_1, lapseTimeMethod)
    # Time step (s) of the Paasschens functions
    DT = 0.5
    lR = length(B0)
    # Kernel grid steps (m) = inversion grid steps / kT (DRz < 0: z descending);
    # DR = smallest step in km, used as distance step of the catalogue
    x1 = sort(unique(modv[:, 1]))
    y1 = sort(unique(modv[:, 2]))
    z1 = sort(unique(modv[:, 3]), rev = true)
    DRx = (x1[2] - x1[1]) / kT
    DRy = (y1[2] - y1[1]) / kT
    DRz = (z1[2] - z1[1]) / kT
    DR = minimum((DRx, DRy, abs(DRz))) / 1000

    # Local reference centred on the source-station midpoint; S = half-size of the
    # kernel grid (m): distance v T travelled at lapse time plus a 5 km cushion
    origin = (event[1:3] .+ station[1:3]) ./ 2
    lapseTimeMethod == "Peak" && (T = mean(T))
    S = (v * T + 5) * 1000
    xyzs = event[1:3] .- origin
    xyzr = station[1:3] .- origin
    # Source-station distance (km)
    D0 = sqrt(sum(abs2, xyzr .- xyzs)) / 1000

    # Grid: horizontally from -S to S; vertically from the top of the inversion
    # grid down to -S
    x_grid = matlab_colon(-S + DRx, DRx, S)
    y_grid = matlab_colon(-S + DRy, DRy, S)
    z_grid = matlab_colon(maximum(z1) - origin[3], DRz, -S)
    Nxyz = length(x_grid) * length(y_grid) * length(z_grid)
    Nxyz > 50e6 && @warn "A lot of points are going to be included in the kernel grid ($Nxyz points)"
    r_grid = unfold_xyz(x_grid, y_grid, z_grid)

    # Distances (km) of each grid point from the source (r1) and the station (r2);
    # only points inside the ellipsoid r1 + r2 < v T get a kernel value
    r1 = sqrt.((r_grid[:, 1] .- xyzs[1]) .^ 2 .+ (r_grid[:, 2] .- xyzs[2]) .^ 2 .+
               (r_grid[:, 3] .- xyzs[3]) .^ 2) ./ 1000
    r2 = sqrt.((r_grid[:, 1] .- xyzr[1]) .^ 2 .+ (r_grid[:, 2] .- xyzr[2]) .^ 2 .+
               (r_grid[:, 3] .- xyzr[3]) .^ 2) ./ 1000
    interior = r1 .+ r2 .< T * v

    # Catalog of Paasschens functions
    # (one per distance R = DR, 2DR, ..., v T)
    R = matlab_colon(DR, DR, T * v)
    NR = length(R)
    t0_R = zeros(NR)
    A_R = zeros(NR, lR)
    codas = Vector{Matrix{Float64}}(undef, NR)
    for (ir, r) in enumerate(R)
        t0, A, _, coda = paasschens(r, v, B0, Le_1, DT, T)
        t0_R[ir] = t0
        A_R[ir, :] .= A
        codas[ir] = coda
    end
    # Pack the codas into NTmax × nFreq × NR (zero-padded). MATLAB preallocates
    # ceil((T+DT)/DT) rows and grows the array if a coda is longer; the max
    # reproduces that
    NTmax = max(ceil(Int, (T + DT) / DT), maximum(size.(codas, 1)))
    codaM = zeros(NTmax, lR, NR)
    for ir in 1:NR
        codaM[1:size(codas[ir], 1), :, ir] .= codas[ir]
    end

    Kint = convolution_catalog(R, t0_R, A_R, codaM, T, DT, v, D0, DR)

    # Kernel on the grid: bilinear interpolation of the catalogue K(r1, r2) at
    # each interior point (exterior points stay 0)
    K = zeros(Nxyz, lR)
    idx = findall(interior)
    for c in 1:lR
        Kc = @view Kint[:, :, c]
        if !isempty(idx)
            for n in idx
                K[n, c] = interp2_linear_clamped(R, R, Kc, r1[n], r2[n])
            end
            # Remaining NaN replaced by the maximum (MATLAB max ignores NaN)
            if any(isnan, view(K, :, c))
                mx = maximum(filter(!isnan, view(K, :, c)))
                K[isnan.(view(K, :, c)), c] .= mx
            end
        end
    end
    # Back to the original reference (m)
    r_out = r_grid .+ reshape(collect(origin), 1, 3)
    return K, r_out
end

"""
    convolution_catalog(R, t0_R, A_R, codaM, T, DT, v, D0, DR) -> Kint

Catalog of the r1 - r2 convolutions of the Paasschens functions evaluated
at lapse time `T` (coda × coda plus coda × delta terms), NR × NR × nFreq.
Entry (i1, i2) is the kernel value for a point at distance `R[i1]` from the source
and `R[i2]` from the station:
`∫ E_coda(R[i1], u) E_coda(R[i2], T - u) du + A(R[i1]) E_coda(R[i2], T - t0(R[i1]))
+ A(R[i2]) E_coda(R[i1], T - t0(R[i2]))`, computed only when
`D0 - 2DR ≤ R[i1] + R[i2] < (T + DT) v` (otherwise 0). Symmetric in (i1, i2).
MATLAB: inlined in `Murat_kernels.m` (`coda_coda + coda_delta`).

# Arguments
- `R`: distances of the catalogue (km), length NR.
- `t0_R`, `A_R`: delta arrival times (s, length NR) and amplitudes (NR × nFreq).
- `codaM`: coda samples, NTmax × nFreq × NR (row u = time `t0 + (u-1) DT`).
- `T`, `DT`: lapse time and time step (s).
- `v`: average velocity (km/s).
- `D0`: source-station distance (km).
- `DR`: distance step (km).

# Returns
- `Kint`: NR × NR × nFreq kernel catalogue.
"""
function convolution_catalog(R, t0_R, A_R, codaM, T, DT, v, D0, DR)
    NR, lR = length(R), size(A_R, 2)
    Kint = zeros(NR, NR, lR)
    # Lower triangle only, mirrored at the end of each step
    for i1 in 1:NR, i2 in 1:i1
        # Number of samples in the overlap of E(r1, t) and E(r2, T - t)
        # (MATLAB round: ties away from zero)
        N_mues = round(Int, (T - t0_R[i1] - t0_R[i2]) / DT, RoundNearestTiesAway)
        if R[i1] + R[i2] >= (D0 - 2DR) && R[i1] + R[i2] < (T + DT) * v && N_mues >= 1
            for c in 1:lR
                # delta1 × coda2 + delta2 × coda1
                cd = A_R[i1, c] * codaM[N_mues+1, c, i2] + A_R[i2, c] * codaM[N_mues+1, c, i1]
                # coda1 × coda2: discrete convolution, sample times add up to ≈ T
                cc = 0.0
                for u in 1:N_mues+1
                    cc += codaM[u, c, i1] * codaM[N_mues+2-u, c, i2]
                end
                Kint[i1, i2, c] = Kint[i2, i1, c] = cc * DT + cd
            end
        end
    end
    return Kint
end

"""
    murat_coda_matrix(modvQc, K, r) -> a

Coda sensitivity kernel moved to the inversion grid: the kernel `K` (one frequency)
is linearly interpolated at the inversion grid nodes `modvQc`; if any value is NaN
or all are zero, NaN and zeros are replaced by `floatmin()` (as MATLAB); the result
is normalized to sum 1. Gives one row of
the Qc inversion matrix (MuRAT list order).
Used with `kernelType = "middle"` (MuRAT 3 behaviour: one kernel, at the middle of
the coda window); the MuRAT 4.0 default is [`murat_coda_matrix_differential`](@ref).
MATLAB: `Murat_codaMatrix.m` of MuRAT 3 (which received two kernels but only used
the second one); the kernel figure is not produced here.

# Arguments
- `modvQc`: inversion grid, nModel × ≥3 `[x y z ...]` in m.
- `K`: kernel values for one frequency (length Nxyz), from `murat_kernels`.
- `r`: Nxyz × 3 kernel grid coordinates in m, from `murat_kernels`.

# Returns
- `a`: Qc inversion matrix row (length nModel), sum 1.
"""
function murat_coda_matrix(modvQc, K::AbstractVector, r)
    # Kernel as meshgrid array (ny, nx, nz); z axis descending as in the list order
    xK = sort(unique(r[:, 1]))
    yK = sort(unique(r[:, 2]))
    zK = sort(unique(r[:, 3]), rev = true)
    KMesh = fold_field(xK, yK, zK, K)
    # Inversion grid nodes as meshgrid arrays
    x = sort(unique(modvQc[:, 1]))
    y = sort(unique(modvQc[:, 2]))
    z = sort(unique(modvQc[:, 3]), rev = true)
    X, Y, Z = meshgrid3(x, y, z)
    # As MATLAB interp3: NaN where the inversion node is outside the kernel grid
    mK = interp3_linear(xK, yK, zK, KMesh, X, Y, Z)
    if any(isnan, mK) || all(==(0), mK)
        mK[isnan.(mK)] .= floatmin(Float64)
        mK[mK .== 0] .= floatmin(Float64)
    end
    # Back to MuRAT list order, normalized
    a = unfold_field(mK)
    return a ./ sum(a)
end

"""
    kernel_on_grid(modvQc, K, r) -> (mK, x, y, z)

Kernel `K` (one frequency, on the kernel grid `r`) linearly interpolated at the
nodes of the inversion grid `modvQc` (NaN outside the kernel grid, as MATLAB
`interp3`).

# Returns
- `mK`: kernel on the inversion grid, meshgrid array (ny, nx, nz).
- `x`, `y`, `z`: inversion grid vectors (m; z descending).
"""
function kernel_on_grid(modvQc, K::AbstractVector, r)
    xK = sort(unique(r[:, 1]))
    yK = sort(unique(r[:, 2]))
    zK = sort(unique(r[:, 3]), rev = true)
    KMesh = fold_field(xK, yK, zK, K)
    x = sort(unique(modvQc[:, 1]))
    y = sort(unique(modvQc[:, 2]))
    z = sort(unique(modvQc[:, 3]), rev = true)
    X, Y, Z = meshgrid3(x, y, z)
    return interp3_linear(xK, yK, zK, KMesh, X, Y, Z), x, y, z
end

"""
    fix_kernel!(mK, K, r, x, y, z) -> mK

Repairs a kernel on the inversion grid, as `fixKernelLocal` in MuRAT 4.0: NaN
values become `floatmin()`, negative values 0; if no value is then above
`floatmin()`, the kernel becomes a single 1 at the inversion node nearest to the
maximum of the original kernel `K` (grid `r`).

# Arguments
- `mK`: kernel on the inversion grid (ny, nx, nz), modified in place.
- `K`, `r`: the original kernel and its grid (from `murat_kernels`).
- `x`, `y`, `z`: inversion grid vectors (m; z descending).
"""
function fix_kernel!(mK, K::AbstractVector, r, x, y, z)
    tiny = floatmin(Float64)
    mK[isnan.(mK)] .= tiny
    mK[mK .< 0] .= 0
    if !any(>(tiny), mK)
        # MATLAB max ignores NaN; ties go to the first index (also for min)
        idx = argmax(map(k -> isnan(k) ? -Inf : k, K))
        nodes = unfold_xyz(x, y, z)
        d = vec(sum((nodes[:, 1:3] .- r[idx:idx, 1:3]) .^ 2, dims = 2))
        single = zeros(size(nodes, 1))
        single[argmin(d)] = 1
        mK .= fold_field(x, y, z, single)
    end
    return mK
end

"""
    murat_coda_matrix_differential(modvQc, KS, rS, KE, rE) -> a

Differential coda sensitivity kernel on the inversion grid (MuRAT 4.0 default,
`kernelType = "differential"`): the kernels at the start (`KS`, lapse time tCoda)
and at the end (`KE`, tCoda + codaWindow) of the coda window are interpolated at
the inversion nodes, each repaired if it has NaN or is all zero, then
`KS - KE` is repaired (`fix_kernel!`: negative values set to 0) and normalized to
sum 1. One row of the Qc inversion matrix (MuRAT list order).
MATLAB: `Murat_codaMatrix.m` of MuRAT 4.0 (`mK = -(mKE - mKS)`); the kernel figure
is not produced here.

# Arguments
- `modvQc`: inversion grid, nModel × ≥3 `[x y z ...]` in m.
- `KS`, `rS`: kernel at the coda start (one frequency) and its grid.
- `KE`, `rE`: kernel at the coda end (one frequency) and its grid.

# Returns
- `a`: Qc inversion matrix row (length nModel), sum 1.
"""
function murat_coda_matrix_differential(modvQc, KS::AbstractVector, rS, KE::AbstractVector, rE)
    mKS, x, y, z = kernel_on_grid(modvQc, KS, rS)
    mKE, _, _, _ = kernel_on_grid(modvQc, KE, rE)
    (any(isnan, mKS) || all(==(0), mKS)) && fix_kernel!(mKS, KS, rS, x, y, z)
    (any(isnan, mKE) || all(==(0), mKE)) && fix_kernel!(mKE, KE, rE, x, y, z)
    # Sensitivity of the energy lost between the start and the end of the window
    mK = -(mKE .- mKS)
    fix_kernel!(mK, KE, rE, x, y, z)
    a = unfold_field(mK)
    return a ./ sum(a)
end
