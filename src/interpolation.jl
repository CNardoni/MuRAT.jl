# Interpolation routines replacing MATLAB interp3, griddedInterpolant,
# griddata (3D linear) and nearest-neighbour NaN filling.

# Locate q in a monotonic (increasing or decreasing) grid vector: returns the
# index i of the interval [g[i], g[i+1]] and the fraction, or 0 if outside.
"""
    grid_locate(g, q) -> (i, f)

Interval of the monotonic (increasing or decreasing) grid vector `g` containing
`q`, and the fractional position of `q` in it. Shared by the linear interpolators
below. MATLAB: no MATLAB counterpart (Julia helper).

# Returns
- `i`: index such that `q` lies between `g[i]` and `g[i+1]` (1 ≤ i ≤ n-1; a point
  on the last node gets i = n-1, f = 1); 0 if `q` is outside `g` or `g` has < 2 nodes.
- `f`: `(q - g[i]) / (g[i+1] - g[i])` in [0, 1]; 0.0 when outside.
"""
function grid_locate(g::AbstractVector, q::Real)
    n = length(g)
    n < 2 && return (0, 0.0)
    if g[end] >= g[1]
        (q < g[1] || q > g[end]) && return (0, 0.0)
        i = clamp(searchsortedlast(g, q), 1, n - 1)
    else
        (q > g[1] || q < g[end]) && return (0, 0.0)
        i = clamp(searchsortedlast(g, q, rev = true), 1, n - 1)
    end
    return (i, (q - g[i]) / (g[i+1] - g[i]))
end

"""
    interp3_linear(x, y, z, V, xq, yq, zq) -> value

MATLAB `interp3(X,Y,Z,V,Xq,Yq,Zq,'linear',NaN)` for a meshgrid field `V`
of size (ny, nx, nz) defined on grid vectors `x`, `y`, `z` (monotonic).
Trilinear interpolation at one query point; NaN outside the grid. Used to move
coda kernels onto the inversion grid (`murat_coda_matrix`).
MATLAB: built-in `interp3` (used in `Murat_codaMatrix.m`).

# Arguments
- `x`, `y`, `z`: grid vectors (length nx, ny, nz), increasing or decreasing.
- `V`: field values, (ny, nx, nz), indexed `V[j, i, k]` as MATLAB meshgrid.
- `xq`, `yq`, `zq`: query point (same units as the grid).

# Returns
- Interpolated value, or NaN if the point is outside the grid.
"""
function interp3_linear(x, y, z, V::AbstractArray{<:Real,3}, xq::Real, yq::Real, zq::Real)
    i, fx = grid_locate(x, xq)
    j, fy = grid_locate(y, yq)
    k, fz = grid_locate(z, zq)
    (i == 0 || j == 0 || k == 0) && return NaN
    # Interpolate along x, then y, then z (V is indexed [y, x, z])
    c00 = V[j, i, k] * (1 - fx) + V[j, i+1, k] * fx
    c10 = V[j+1, i, k] * (1 - fx) + V[j+1, i+1, k] * fx
    c01 = V[j, i, k+1] * (1 - fx) + V[j, i+1, k+1] * fx
    c11 = V[j+1, i, k+1] * (1 - fx) + V[j+1, i+1, k+1] * fx
    c0 = c00 * (1 - fy) + c10 * fy
    c1 = c01 * (1 - fy) + c11 * fy
    return c0 * (1 - fz) + c1 * fz
end

"""
    interp3_linear(x, y, z, V, Xq, Yq, Zq) -> array of size(Xq)

Array version of `interp3_linear`: evaluates the trilinear interpolant at every
query point of the arrays `Xq`, `Yq`, `Zq` (same size, e.g. from `meshgrid3`) and
returns an array of the same shape (NaN outside the grid).
"""
interp3_linear(x, y, z, V, Xq::AbstractArray, Yq::AbstractArray, Zq::AbstractArray) =
    [interp3_linear(x, y, z, V, Xq[n], Yq[n], Zq[n]) for n in eachindex(Xq)] |> v -> reshape(v, size(Xq))

"""
    interp2_linear_clamped(r1, r2, K, q1, q2) -> value

MATLAB `griddedInterpolant({r1,r2}, K, 'linear', 'nearest')` evaluated at
(q1, q2): bilinear inside the grid; outside it (in any dimension), the
value of the nearest grid node.
Used to read the kernel catalogue K(r1, r2) at the source/station distances of
each kernel grid point (`murat_kernels`).
MATLAB: built-in `griddedInterpolant` (used in `Murat_kernels.m`).

# Arguments
- `r1`, `r2`: increasing grid vectors (here distances in km).
- `K`: values, length(r1) × length(r2).
- `q1`, `q2`: query point.

# Returns
- Interpolated (or nearest-node) value.
"""
function interp2_linear_clamped(r1::AbstractVector, r2::AbstractVector, K::AbstractMatrix, q1::Real, q2::Real)
    if q1 < r1[1] || q1 > r1[end] || q2 < r2[1] || q2 > r2[end]
        return K[argmin(abs.(r1 .- q1)), argmin(abs.(r2 .- q2))]
    end
    i, f1 = grid_locate(r1, q1)
    j, f2 = grid_locate(r2, q2)
    return (K[i, j] * (1 - f1) + K[i+1, j] * f1) * (1 - f2) +
           (K[i, j+1] * (1 - f1) + K[i+1, j+1] * f1) * f2
end

"""
    griddata3_linear(xo, yo, zo, vo, Xq, Yq, Zq) -> array of size(Xq)

MATLAB `griddata(xo,yo,zo,vo,Xq,Yq,Zq)` (linear, 3D): barycentric
interpolation on the Delaunay tetrahedralization of the scattered points;
NaN outside the convex hull.
Used by `rescale_field` to resample an input 3D velocity model onto the MuRAT grids.
MATLAB: built-in `griddata` (used in `Murat_rescale.m`). A query point lying on a
face shared by two tetrahedra takes the value from the first one found; the value
is the same up to rounding.

# Arguments
- `xo`, `yo`, `zo`, `vo`: scattered point coordinates and values (vectors).
- `Xq`, `Yq`, `Zq`: query coordinates (arrays of equal size).

# Returns
- Interpolated values with the shape of `Xq`; NaN outside the convex hull.
"""
function griddata3_linear(xo, yo, zo, vo, Xq::AbstractArray, Yq::AbstractArray, Zq::AbstractArray)
    # Points as a 3 × nPoints matrix; tets is 4 × nTet (vertex indices)
    P = permutedims(hcat(float.(xo), float.(yo), float.(zo)))
    tets = MiniQhull.delaunay(P)        # qhull d Qt Qbb Qc Qz
    ntet = size(tets, 2)

    # Bin tetrahedra on a coarse grid of their bounding boxes
    # (nb³ bins over the data bounding box, ~4 tetrahedra per bin on average;
    # each tetrahedron is listed in every bin its bounding box touches)
    lo = vec(minimum(P, dims = 2))
    hi = vec(maximum(P, dims = 2))
    nb = max(1, round(Int, cbrt(ntet / 4)))
    # Bin index (1..nb) of coordinate v along dimension d
    cell(v, d) = clamp(1 + floor(Int, (v - lo[d]) / max(hi[d] - lo[d], eps()) * nb), 1, nb)
    bins = [Int[] for _ in 1:nb, _ in 1:nb, _ in 1:nb]
    for t in 1:ntet
        idx = tets[:, t]
        bmin = [minimum(P[d, idx]) for d in 1:3]
        bmax = [maximum(P[d, idx]) for d in 1:3]
        for a in cell(bmin[1], 1):cell(bmax[1], 1), b in cell(bmin[2], 2):cell(bmax[2], 2),
            c in cell(bmin[3], 3):cell(bmax[3], 3)
            push!(bins[a, b, c], t)
        end
    end

    out = fill(NaN, size(Xq))
    for n in eachindex(Xq)
        q = (Xq[n], Yq[n], Zq[n])
        # Quick reject: outside the data bounding box (small relative tolerance)
        any(q[d] < lo[d] - 1e-9 * (hi[d] - lo[d]) || q[d] > hi[d] + 1e-9 * (hi[d] - lo[d]) for d in 1:3) && continue
        for t in bins[cell(q[1], 1), cell(q[2], 2), cell(q[3], 3)]
            # Barycentric coordinates: solve M λ = q - p1 with M = edges from p1;
            # skip degenerate (flat) tetrahedra
            v = tets[:, t]
            p1 = P[:, v[1]]
            M = hcat(P[:, v[2]] .- p1, P[:, v[3]] .- p1, P[:, v[4]] .- p1)
            detM = det(M)
            abs(detM) < 1e-12 * (maximum(abs.(M)) + 1)^3 && continue
            λ = M \ ([q...] .- p1)
            l1 = 1 - sum(λ)
            # Inside (or on the boundary of) the tetrahedron if all weights ≥ 0
            tolb = -1e-10
            if all(λ .>= tolb) && l1 >= tolb
                out[n] = l1 * vo[v[1]] + λ[1] * vo[v[2]] + λ[2] * vo[v[3]] + λ[3] * vo[v[4]]
                break
            end
        end
    end
    return out
end

"""
    rescale_field(xo, yo, zo, vo, x, y, z) -> (X, Y, Z, list, mesh)

`Murat_rescale`: scattered field to the meshgrid of `x`, `y`, `z` (linear,
NaN gaps filled with `fill_nan3`). Returns `(X, Y, Z, list, mesh)`.
Used by `murat_modv3d` to put an input 3D velocity model on the inversion and
propagation grids.
MATLAB: `Murat_rescale.m`.

# Arguments
- `xo`, `yo`, `zo`, `vo`: scattered points (m) and values (e.g. velocity, km/s).
- `x`, `y`, `z`: target grid vectors (m).

# Returns
- `X`, `Y`, `Z`: meshgrid arrays of the target grid, (ny, nx, nz).
- `list`: nx*ny*nz × 4 `[x y z v]` in MuRAT list order (`unfold_mesh`).
- `mesh`: interpolated values as a (ny, nx, nz) array.
"""
function rescale_field(xo, yo, zo, vo, x, y, z)
    X, Y, Z = meshgrid3(x, y, z)
    mR = griddata3_linear(xo, yo, zo, vo, X, Y, Z)
    if any(isnan, mR)
        mR = fill_nan3(mR)
    end
    return X, Y, Z, unfold_mesh(X, Y, Z, mR), mR
end

"""
    fill_nearest!(v) -> v

Fills the NaN of a vector with the nearest non-NaN value, as MATLAB
`fillmissing(v, 'nearest')`: at equal distance the later value is used; a vector
with no valid value stays NaN.
"""
function fill_nearest!(v::AbstractVector)
    valid = findall(!isnan, v)
    isempty(valid) && return v
    for i in eachindex(v)
        isnan(v[i]) || continue
        k = searchsortedfirst(valid, i)            # first valid index after i
        if k > length(valid)
            v[i] = v[valid[end]]
        elseif k == 1
            v[i] = v[valid[1]]
        else
            before, after = valid[k-1], valid[k]
            v[i] = (i - before) < (after - i) ? v[before] : v[after]
        end
    end
    return v
end

"""
    fill_nan3(A) -> array

Fills the NaN of a 3D array by nearest-neighbour propagation along each dimension
(two passes over dimensions 1, 2, 3), as `fillNaN3` in `Murat_rescale.m` of
MuRAT 4.0 (`fillmissing(A, 'nearest', dim)`). Used where the velocity model is
interpolated outside the convex hull of its nodes.
MATLAB: `fillNaN3` in `Murat_rescale.m` (MuRAT 3 used `inpaintn`).
"""
function fill_nan3(A::AbstractArray{<:Real,3})
    A = float.(A)
    for pass in 1:2
        for dim in 1:3
            for line in eachslice(A; dims = Tuple(d for d in 1:3 if d != dim))
                fill_nearest!(line)
            end
        end
        any(isnan, A) || return A
    end
    return A
end
