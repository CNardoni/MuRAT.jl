# General helpers reproducing the MATLAB behaviour MuRAT relies on.

"""
    matlab_colon(a, d, b) -> Vector{Float64}

Vector `a:d:b` computed exactly as MATLAB's colon operator (same number of
elements and same values, built symmetrically from both ends).
Used wherever MuRAT builds ranges with `:` (kernel grids, Paasschens time vectors,
cell boundaries) so that lengths and values match MATLAB.
MATLAB: built-in colon operator `a:d:b`. Empty if `d == 0` or the step points away
from `b`.

Caveat: the tolerance used to decide whether `b` is reached is
`2 eps max(|a|, |b|)`, added to the element count `(b - a)/d` without dividing by
`d`. When the step is small relative to the endpoints (e.g. `t + si : si : t + n*si`
with a sampling interval `si` well below 1), the rounding error of `(b - a)/d`
(about `eps |b| / d`) can exceed it and the last element is dropped (n - 1 elements
instead of n). When an exact count is
needed, build the range from an integer count instead (as done in `murat_qc`).

# Arguments
- `a`, `d`, `b`: start, step and end (converted to Float64).

# Returns
- The vector of values from `a` to (at most) `b` with step `d`.
"""
function matlab_colon(a::Real, d::Real, b::Real)
    a, d, b = float(a), float(d), float(b)
    (d == 0 || (d > 0 && a > b) || (d < 0 && a < b)) && return Float64[]
    tol = 2.0 * eps() * max(abs(a), abs(b))
    sig = sign(d)
    # n = number of intervals (length - 1); integer start with unit step is exact
    n = (a == floor(a) && d == 1) ? Int(floor(b) - a) : Int(floor((b - a) / d + tol))
    # Right end point, snapped to b if within tolerance
    c = a + n * d
    if sig * (c - b) > -tol
        c = b
    end
    # First half from a forwards, second half from c backwards (MATLAB symmetry)
    v = zeros(n + 1)
    for k in 0:fld(n, 2)
        v[1+k] = a + k * d
        v[n+1-k] = c - k * d
    end
    if iseven(n)
        v[n÷2+1] = (a + c) / 2
    end
    return v
end

"""
    smooth_moving(y, span) -> Vector{Float64}

MATLAB `smooth(y, span)` with the default moving average: the span is
truncated to an odd integer and the window shrinks symmetrically at the ends.
Used by MuRAT 3 to smooth the squared Hilbert envelopes; MuRAT 4.0 uses RMS
envelopes (`rms_envelope`), so it is no longer called by the package.
MATLAB: `smooth` (Curve Fitting Toolbox), called in `Murat_envelope.m` of MuRAT 3.

# Arguments
- `y`: signal samples.
- `span`: window length in samples (non-integer allowed; capped at length(y)).

# Returns
- Smoothed signal, same length as `y`.
"""
function smooth_moving(y::AbstractVector, span::Real)
    n = length(y)
    s = min(floor(Int, span), n)
    # Largest odd width ≤ s
    width = s - 1 + mod(s, 2)
    width <= 1 && return float.(collect(y))
    out = similar(y, Float64)
    half = (width - 1) ÷ 2
    # Moving mean via cumulative sum; near the ends the half-width k shrinks so
    # that the window stays centred
    cs = cumsum(vcat(0.0, float.(y)))
    for i in 1:n
        k = min(half, i - 1, n - i)
        out[i] = (cs[i+k+1] - cs[i-k]) / (2k + 1)
    end
    return out
end

"""
    detrend_linear(x) -> Vector{Float64}

Remove the least-squares straight line (MATLAB `detrend(x,1)`).
Applied to each waveform before tapering and filtering (`murat_envelope`).
MATLAB: `detrend(x, 1)`, called in `Murat_envelope.m`.
"""
function detrend_linear(x::AbstractVector)
    n = length(x)
    t = collect(1.0:n)
    G = hcat(ones(n), t)
    c = G \ float.(x)
    return float.(x) .- G * c
end

"""
    tukeywin(n, r) -> Vector{Float64}

Tukey (tapered cosine) window, MATLAB `tukeywin(n, r)`.
`r` is the fraction of the window inside the cosine tapers: `r ≤ 0` gives a
rectangular window, `r ≥ 1` a Hann window. Used as a 5% taper on each waveform
(`murat_envelope`).
MATLAB: `tukeywin` (Signal Processing Toolbox), called in `Murat_envelope.m`.
"""
function tukeywin(n::Integer, r::Real)
    r <= 0 && return ones(n)
    r >= 1 && return 0.5 .* (1 .- cos.(2π .* (0:n-1) ./ (n - 1)))
    t = range(0, 1, length = n)
    per = r / 2
    # tl = last sample of the left taper, th = first sample of the right taper
    tl = floor(Int, per * (n - 1)) + 1
    th = n - tl + 1
    w = ones(n)
    for i in 1:tl
        w[i] = 0.5 * (1 + cos(π / per * (t[i] - per)))
    end
    for i in th:n
        w[i] = 0.5 * (1 + cos(π / per * (t[i] - 1 + per)))
    end
    return w
end

"""
    trapz_columns(Y, dt) -> Vector{Float64}

Trapezoidal integral of each column of `Y` with uniform spacing `dt`.
Used to integrate envelope energy in the body, coda and noise windows
(`murat_body`). Returns zeros if `Y` has fewer than 2 rows.
MATLAB: `trapz(0:dt:..., Y)`, called in `Murat_body.m`.
"""
function trapz_columns(Y::AbstractMatrix, dt::Real)
    size(Y, 1) < 2 && return zeros(size(Y, 2))
    # dt * (sum of all samples - half of the first and last)
    return vec(dt .* (sum(Y, dims = 1) .- (Y[1:1, :] .+ Y[end:end, :]) ./ 2))
end

# --- MuRAT grid formats ---------------------------------------------------
# MuRAT stores 3D fields as lists ordered with x slowest and z fastest
# ("unfolded"), or as meshgrid arrays of size (ny, nx, nz) ("folded").

"""
    unfold_xyz(x, y, z) -> Matrix (nx*ny*nz × 3)

Unfold grid vectors into MuRAT list format `[x y z]` (x slowest, z fastest).
Used to build the inversion grid (`checks.jl`), the kernel grid (`murat_kernels`)
and the test-input grids (`inversion.jl`).
MATLAB: `Murat_unfoldXYZ.m`.
"""
function unfold_xyz(x, y, z)
    lx, ly, lz = length(x), length(y), length(z)
    r = zeros(lx * ly * lz, 3)
    idx = 0
    for i in 1:lx, j in 1:ly, k in 1:lz
        idx += 1
        r[idx, 1] = x[i]
        r[idx, 2] = y[j]
        r[idx, 3] = z[k]
    end
    return r
end

"""
    meshgrid3(x, y, z) -> (X, Y, Z)

MATLAB `meshgrid(x, y, z)`: arrays of size (ny, nx, nz).
`X[j, i, k] = x[i]`, `Y[j, i, k] = y[j]`, `Z[j, i, k] = z[k]`.
MATLAB: built-in `meshgrid` (also the X, Y, Z outputs of `Murat_fold.m`).
"""
function meshgrid3(x, y, z)
    lx, ly, lz = length(x), length(y), length(z)
    X = [x[i] for j in 1:ly, i in 1:lx, k in 1:lz]
    Y = [y[j] for j in 1:ly, i in 1:lx, k in 1:lz]
    Z = [z[k] for j in 1:ly, i in 1:lx, k in 1:lz]
    return X, Y, Z
end

"""
    fold_field(x, y, z, v) -> Array (ny, nx, nz)

`Murat_fold`: field `v` in MuRAT list order to a meshgrid array (ny, nx, nz).
Only the lengths of `x`, `y`, `z` are used. Unlike MATLAB, the X, Y, Z meshgrid
arrays are not returned (use `meshgrid3`).
MATLAB: `Murat_fold.m`.

# Arguments
- `x`, `y`, `z`: grid vectors (lengths nx, ny, nz).
- `v`: values in MuRAT list order (x slowest, z fastest), length nx*ny*nz.
"""
function fold_field(x, y, z, v)
    lx, ly, lz = length(x), length(y), length(z)
    # List order -> (nz, ny, nx) column-major array -> permute to (ny, nx, nz)
    return permutedims(reshape(collect(v), lz, ly, lx), (2, 3, 1))
end

"""
    unfold_field(V) -> Vector

Inverse of `fold_field`: meshgrid array (ny, nx, nz) to MuRAT list order.
MATLAB: `reshape(permute(V, [3 1 2]), [], 1)` (e.g. in `Murat_codaMatrix.m`).
"""
unfold_field(V::AbstractArray{<:Any,3}) = vec(permutedims(V, (3, 1, 2)))

"""
    unfold_mesh(X, Y, Z, V) -> Matrix (nx*ny*nz × 4)

`Murat_unfold`: meshgrid arrays to an n×4 list `[x y z v]` in MuRAT order.
All inputs are (ny, nx, nz) arrays, e.g. from `meshgrid3`.
MATLAB: `Murat_unfold.m`.
"""
function unfold_mesh(X, Y, Z, V)
    ly, lx, lz = size(X)
    out = zeros(lx * ly * lz, 4)
    idx = 0
    for i in 1:lx, j in 1:ly, k in 1:lz
        idx += 1
        out[idx, :] .= (X[j, i, k], Y[j, i, k], Z[j, i, k], V[j, i, k])
    end
    return out
end

"""
    dd_coordinates(origin, ending, nLat, nLong, nzc) -> Matrix (nLat*nLong*nzc × 3)

`Murat_DDcoordinates`: grid nodes in degrees, x (longitude) fastest.
Rows are `[longitude latitude depth]`, with depth slowest. Note this order differs
from the MuRAT list order used elsewhere (x slowest, z fastest).
MATLAB: `Murat_DDcoordinates.m`.

# Arguments
- `origin`, `ending`: grid corners `[lat, lon, z]` (degrees, degrees, z units of the
  settings).
- `nLat`, `nLong`, `nzc`: number of nodes in latitude, longitude and depth.
"""
function dd_coordinates(origin, ending, nLat, nLong, nzc)
    yLat = range(origin[1], ending[1], length = nLat)
    xLong = range(origin[2], ending[2], length = nLong)
    zDep = range(origin[3], ending[3], length = nzc)
    out = zeros(nLat * nLong * nzc, 3)
    idx = 0
    # Same order as MATLAB ndgrid(xLong, yLat, zDep)(:)
    for k in 1:nzc, j in 1:nLat, i in 1:nLong
        idx += 1
        out[idx, :] .= (xLong[i], yLat[j], zDep[k])
    end
    return out
end

# --- MATLAB-style sorting and uniqueness ------------------------------------

"""
    sortrows_matlab(A, cols) -> (sorted A, permutation)

MATLAB `sortrows(A, cols)`: stable sort by the given columns, a negative
column number sorting in descending order.
Used in `murat_declustering` and `murat_modv3d`.
MATLAB: built-in `sortrows` (e.g. in `Murat_declustering.m`, `Murat_modv3D.m`).
"""
function sortrows_matlab(A::AbstractMatrix, cols::AbstractVector{<:Integer})
    p = collect(1:size(A, 1))
    # Successive stable sorts from the last key to the first = lexicographic sort
    for c in reverse(cols)
        col = abs(c)
        sort!(p, by = i -> A[i, col], rev = c < 0, alg = Base.Sort.MergeSort)
    end
    return A[p, :], p
end

"""
    unique_rows_stable(A) -> Matrix

Unique rows of `A` in order of first appearance (MATLAB `'rows','stable'`).
MATLAB: `unique(A, 'rows', 'stable')`, in `Murat_declustering.m`.
"""
function unique_rows_stable(A::AbstractMatrix)
    seen = Set{Vector{Float64}}()
    keep = Int[]
    for i in axes(A, 1)
        r = Vector{Float64}(A[i, :])
        if !(r in seen)
            push!(seen, r)
            push!(keep, i)
        end
    end
    return A[keep, :]
end

"""
    unique_last_indices(v) -> Vector{Int}

Indices of the last occurrence of each unique value, sorted by value
(MATLAB `unique(v,'last')`).
MATLAB: second output of `unique(v, 'rows', 'last')`, in `Murat_declustering.m`.
"""
function unique_last_indices(v::AbstractVector)
    last = Dict{eltype(v),Int}()
    # Later occurrences overwrite earlier ones
    for (i, x) in pairs(v)
        last[x] = i
    end
    vals = sort(collect(keys(last)))
    return [last[x] for x in vals]
end

# --- Synthetic test inputs ----------------------------------------------------

"""
    checkerboard3d(siz, blockSize = 1) -> BitArray

3D checkerboard of size `siz` with cubic blocks of `blockSize` cells: a cell is
`true` when the sum of its block indices `ceil(i/blockSize)` (in the three
dimensions) is odd. Gives the input pattern of the checkerboard tests.
MATLAB: `Murat_checkerboard3D.m` of MuRAT 4.0 (MIT replacement of GIBBON
`checkerBoard3D`, which gives the same pattern).

# Arguments
- `siz`: `(n1, n2, n3)` (missing dimensions are 1).
- `blockSize`: block edge length in cells.
"""
function checkerboard3d(siz, blockSize::Integer = 1)
    n = (collect(siz)..., 1, 1)[1:3]
    return BitArray(isodd(cld(i, blockSize) + cld(j, blockSize) + cld(k, blockSize))
                    for i in 1:n[1], j in 1:n[2], k in 1:n[3])
end

# --- Robust straight-line fit -------------------------------------------------

"""
    robust_line_fit(x, y; maxiter = 50, D = 1e-6) -> [p1, p2]

Straight line `y = p1*x + p2` with bisquare robust weights, following
MATLAB `fit(x, y, 'poly1', 'Robust', 'on')` (iRobustLinearFit): residuals
adjusted by the leverage, scale from the median of the absolute adjusted
residuals without the smallest p-1, tuning constant 4.685. Returns `[p1, p2]`.
Used to fit log10(peak delay) against travel time in the peak-delay selection
(`murat_retain_peak_delay`).
MATLAB: `fit(..., 'poly1', 'Robust', 'on')` (Curve Fitting Toolbox), in
`Murat_retainPeakDelay.m`.

# Arguments
- `x`, `y`: data vectors.
- `maxiter`: maximum number of reweighting iterations.
- `D`: relative convergence tolerance on the coefficients.
"""
function robust_line_fit(x::AbstractVector, y::AbstractVector; maxiter = 50, D = 1e-6)
    y = float.(collect(y))
    n = length(x)
    G = hcat(float.(x), ones(n))
    np = 2
    # Leverages h from the thin QR of G; residual adjustment 1/sqrt(1 - h)
    Q = Matrix(qr(G).Q)[:, 1:np]
    h = min.(0.9999, vec(sum(Q .^ 2, dims = 2)))
    adj = 1 ./ sqrt.(1 .- h)
    tiny_s = 1e-6 * std(y)
    tiny_s == 0 && (tiny_s = 1.0)
    # Start from the ordinary least-squares fit
    p = G \ y
    res = G * p .- y
    iter = 0
    while true
        iter += 1
        # Robust scale S = MAD/0.6745 of the adjusted residuals, bisquare weights
        radj = res .* adj
        sar = sort(abs.(radj))
        S = max(tiny_s, median(sar[np:end]) / 0.6745)
        u = radj ./ (S * 4.685)
        w = [abs(ui) <= 1 ? (1 - ui^2)^2 : 0.0 for ui in u]
        all(==(0), w) && (w .= 1)
        # Weighted least squares
        p0 = p
        sw = sqrt.(w)
        p = (G .* sw) \ (y .* sw)
        res = G * p .- y
        all(abs.(p .- p0) .<= D .* max.(abs.(p), abs.(p0))) && break
        iter > maxiter && break
    end
    return [p[1], p[2]]
end

"""
    rms_envelope(x, n) -> Vector

Upper RMS envelope of `x` over a sliding window of `n` samples, as MATLAB
`envelope(x, n, 'rms')`: the mean of `x` is removed, the RMS is computed over the
window (for even `n`, n/2 - 1 samples before and n/2 after each sample; the
window shrinks at the ends) and the mean is added back.
MATLAB: `envelope` (Signal Processing Toolbox), used in `Murat_envelope.m`.
"""
function rms_envelope(x::AbstractVector, n::Integer)
    xm = sum(x) / length(x)
    c = [0.0; cumsum(abs2.(x .- xm))]
    nx = length(x)
    before = iseven(n) ? n ÷ 2 - 1 : n ÷ 2
    after = n ÷ 2
    y = similar(x, Float64)
    for i in 1:nx
        k1, k2 = max(i - before, 1), min(i + after, nx)
        y[i] = sqrt(max(c[k2+1] - c[k1], 0.0) / (k2 - k1 + 1))
    end
    return xm .+ y
end

"""
    pchip(x, y, xq) -> Vector

Shape-preserving piecewise cubic Hermite interpolation of the points `(x, y)`
(`x` strictly increasing) at `xq`, as MATLAB `interp1(x, y, xq, 'pchip', 'extrap')`:
interior slopes are weighted harmonic means of the neighbouring secant slopes (0
where they change sign); end slopes use the three-point formula of MATLAB's
`pchip`, limited to keep the shape; points outside `x` use the end pieces.
MATLAB: `interp1(..., 'pchip', 'extrap')` (used in `Murat_readLithos1.m`).
"""
function pchip(x::AbstractVector, y::AbstractVector, xq)
    n = length(x)
    T = promote_type(eltype(y), Float32)
    h = diff(x)
    del = diff(y) ./ h
    d = zeros(T, n)
    if n == 2
        d .= del[1]
    else
        # Interior slopes (Fritsch-Carlson, weighted harmonic mean)
        for k in 2:n-1
            if del[k-1] * del[k] > 0
                w1, w2 = 2h[k] + h[k-1], h[k] + 2h[k-1]
                d[k] = (w1 + w2) / (w1 / del[k-1] + w2 / del[k])
            end
        end
        # End slopes: non-centred three-point formula, shape-preserving limits
        function endslope(h1, h2, d1, d2)
            s = ((2h1 + h2) * d1 - h1 * d2) / (h1 + h2)
            if sign(s) != sign(d1)
                s = zero(s)
            elseif sign(d1) != sign(d2) && abs(s) > abs(3d1)
                s = 3d1
            end
            return s
        end
        d[1] = endslope(h[1], h[2], del[1], del[2])
        d[n] = endslope(h[n-1], h[n-2], del[n-1], del[n-2])
    end
    out = similar(xq, T)
    for (i, q) in enumerate(xq)
        k = clamp(searchsortedlast(x, q), 1, n - 1)
        s = q - x[k]
        # Hermite cubic on interval k in powers of s (as MATLAB pwch)
        c = (3del[k] - 2d[k] - d[k+1]) / h[k]
        b = (d[k] - 2del[k] + d[k+1]) / h[k]^2
        out[i] = y[k] + s * (d[k] + s * (c + s * b))
    end
    return out
end
