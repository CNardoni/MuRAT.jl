# Pseudo-bending ray tracing and ray-segment inversion matrices
# (Murat_rays, Murat_tracing, Murat_bending, Murat_traveling, Murat_velocity,
# Murat_veldicity, Murat_cornering, Murat_segments).
#
# Conventions shared by this file:
# - `grid` is the propagation grid `(x = ..., y = ..., z = ...)` built in `checks.jl`:
#   node coordinates in m from the grid origin; x and y increasing, z = altitude in m
#   (negative below sea level) sorted in DECREASING order (shallowest node first).
# - `pvel` is the velocity on that grid as a meshgrid array of size (ny, nx, nz),
#   indexed `pvel[j, i, k]` (y index first, as MATLAB), in km/s.
# - A cell is identified by its shallowest south-west node (ip, jp, kp): it spans
#   x[ip]..x[ip+1], y[jp]..y[jp+1], z[kp] (top)..z[kp+1] (bottom).

# Cell of the propagation grid containing a point; flag = 1 if outside.
"""
    murat_cornering(xx, yy, zz, grid) -> (ip, jp, kp, flag)

Finds the cell of the propagation grid that contains the point (`xx`, `yy`, `zz`),
i.e. the indices of its shallowest south-west corner node. Used by the velocity
interpolation during ray tracing.
MATLAB: `Murat_cornering.m`. Difference: if the point lies exactly on the first node
of an axis, MATLAB `find` would return empty; here the index falls back to 1.

# Arguments
- `xx`, `yy`, `zz`: point coordinates in m (z = altitude, negative downwards).
- `grid`: propagation grid with fields `x`, `y`, `z` (m), z in decreasing order.

# Returns
- `ip`, `jp`, `kp`: indices into `grid.x`, `grid.y`, `grid.z` of the cell corner
  (smallest x, smallest y, largest z).
- `flag`: 0 inside the grid; 1 if the point is outside, in which case the indices are
  those of the nearest node (clamped to the violated boundary), not of a cell.
"""
function murat_cornering(xx, yy, zz, grid)
    xG, yG, zG = grid.x, grid.y, grid.z
    # Outside the grid: checked in the MATLAB order z, then y, then x; return the
    # nearest node along the other two axes and the boundary node along this one.
    if zz < minimum(zG)
        return argmin(abs.(xx .- xG)), argmin(abs.(yy .- yG)), argmin(zG), 1
    elseif zz > maximum(zG)
        return argmin(abs.(xx .- xG)), argmin(abs.(yy .- yG)), argmax(zG), 1
    end
    if yy < minimum(yG)
        return argmin(abs.(xx .- xG)), argmin(yG), argmin(abs.(zz .- zG)), 1
    elseif yy > maximum(yG)
        return argmin(abs.(xx .- xG)), argmax(yG), argmin(abs.(zz .- zG)), 1
    end
    if xx < minimum(xG)
        return argmin(xG), argmin(abs.(yy .- yG)), argmin(abs.(zz .- zG)), 1
    elseif xx > maximum(xG)
        return argmax(xG), argmin(abs.(yy .- yG)), argmin(abs.(zz .- zG)), 1
    end
    # Inside: last node west/south of the point and, since z decreases, the last
    # node above it (shallowest-SW-corner convention of MuRAT).
    ip = something(findlast(<(xx), xG), 1)
    jp = something(findlast(<(yy), yG), 1)
    kp = something(findlast(>(zz), zG), 1)
    return ip, jp, kp, 0
end

# Trilinear interpolation weights inside the cell of a point
"""
    cell_fractions(xx, yy, zz, grid, ip, jp, kp) -> (xf, yf, zf, xd, yd, zd)

Fractional position of a point inside the grid cell with corner (`ip`, `jp`, `kp`),
used as trilinear interpolation weights by `murat_velocity` and `murat_veldicity`.
MATLAB: no separate counterpart (inlined in `Murat_velocity.m`/`Murat_veldicity.m`).

# Returns
- `xf`, `yf`, `zf`: fractions in [0, 1] from the corner node along x, y, z.
- `xd`, `yd`, `zd`: cell sizes in m (`zd` is negative, as z decreases with index).
"""
function cell_fractions(xx, yy, zz, grid, ip, jp, kp)
    xG, yG, zG = grid.x, grid.y, grid.z
    xd = xG[ip+1] - xG[ip]
    yd = yG[jp+1] - yG[jp]
    zd = zG[kp+1] - zG[kp]
    return (xx - xG[ip]) / xd, (yy - yG[jp]) / yd, (zz - zG[kp]) / zd, xd, yd, zd
end

"""
    murat_velocity(xx, yy, zz, grid, pvel) -> v

Velocity at a point by trilinear interpolation of the propagation model `pvel`
(size (ny, nx, nz)). Outside the grid it returns the value of the nearest node.
MATLAB: `Murat_velocity.m`.

# Arguments
- `xx`, `yy`, `zz`: point coordinates in m.
- `grid`: propagation grid (see `murat_cornering`).
- `pvel`: velocity on the grid, (ny, nx, nz); km/s as stored, but callers sometimes
  pass it in m/s (see `murat_traveling`): the result has the units of `pvel`.

# Returns
- Interpolated velocity, same units as `pvel`.
"""
function murat_velocity(xx, yy, zz, grid, pvel)
    ip, jp, kp, flag = murat_cornering(xx, yy, zz, grid)
    flag > 0 && return pvel[jp, ip, kp]
    xf, yf, zf, _, _, _ = cell_fractions(xx, yy, zz, grid, ip, jp, kp)
    # Interpolate along x on the four cell edges, then along y, then along z
    c00 = pvel[jp, ip, kp] * (1 - xf) + pvel[jp, ip+1, kp] * xf
    c10 = pvel[jp+1, ip, kp] * (1 - xf) + pvel[jp+1, ip+1, kp] * xf
    c01 = pvel[jp, ip, kp+1] * (1 - xf) + pvel[jp, ip+1, kp+1] * xf
    c11 = pvel[jp+1, ip, kp+1] * (1 - xf) + pvel[jp+1, ip+1, kp+1] * xf
    c0 = c00 * (1 - yf) + c10 * yf
    c1 = c01 * (1 - yf) + c11 * yf
    return c0 * (1 - zf) + c1 * zf
end

"""
    murat_veldicity(xx, yy, zz, grid, pvel) -> (vx, vy, vz)

Velocity gradient at a point: analytic derivative of the trilinear interpolant
used by `murat_velocity`. Needed by the bending step to know in which direction
the ray should be pushed (towards higher velocity). Zero gradient outside the grid.
MATLAB: `Murat_veldicity.m`.

# Arguments
- `xx`, `yy`, `zz`: point coordinates in m.
- `grid`, `pvel`: as in `murat_velocity`.

# Returns
- `vx`, `vy`, `vz`: dv/dx, dv/dy, dv/dz in (units of `pvel`) per m.
"""
function murat_veldicity(xx, yy, zz, grid, pvel)
    ip, jp, kp, flag = murat_cornering(xx, yy, zz, grid)
    flag > 0 && return 0.0, 0.0, 0.0
    xf, yf, zf, xd, yd, zd = cell_fractions(xx, yy, zz, grid, ip, jp, kp)
    xf1, yf1, zf1 = 1 - xf, 1 - yf, 1 - zf
    # Corner values: vXYZ with 0/1 = lower/upper index along x, y, z
    v000 = pvel[jp, ip, kp];     v100 = pvel[jp, ip+1, kp]
    v010 = pvel[jp+1, ip, kp];   v110 = pvel[jp+1, ip+1, kp]
    v001 = pvel[jp, ip, kp+1];   v101 = pvel[jp, ip+1, kp+1]
    v011 = pvel[jp+1, ip, kp+1]; v111 = pvel[jp+1, ip+1, kp+1]
    # Each component: corner differences along one axis, weighted by the position
    # along the other two, divided by the cell size along that axis
    vx = (yf1 * zf1 * (v100 - v000) + yf * zf1 * (v110 - v010) +
          yf1 * zf * (v101 - v001) + yf * zf * (v111 - v011)) / xd
    vy = (xf1 * zf1 * (v010 - v000) + xf * zf1 * (v110 - v100) +
          xf1 * zf * (v011 - v001) + xf * zf * (v111 - v101)) / yd
    vz = (xf1 * yf1 * (v001 - v000) + xf * yf1 * (v101 - v100) +
          xf1 * yf * (v011 - v010) + xf * yf * (v111 - v110)) / zd
    return vx, vy, vz
end

"""
    murat_traveling(x, y, z, grid, v, pvel) -> (ta, tra)

Cumulative travel time along the ray polyline (`x`, `y`, `z`), used by
`murat_tracing` to decide when the pseudo-bending has converged.
MATLAB: `Murat_traveling.m` (reproduced as is, including its unit mix, see below).

Two branches, chosen from the length of the FIRST segment only:
- short segments (≤ 5 propagation cells): time from the velocities `v` at the ray
  points (km/s) and lengths in m, so the time is in ms;
- long segments: each segment is split into pieces ≤ 5 cells and the velocity is
  sampled at their midpoints from `pvel` converted to m/s, so the time is in s.

# Arguments
- `x`, `y`, `z`: ray point coordinates in m (vectors of length n ≥ 2).
- `grid`: propagation grid (see `murat_cornering`).
- `v`: velocity at each ray point (km/s), from `murat_velocity`.
- `pvel`: propagation velocity model (ny, nx, nz), km/s.

# Returns
- `ta`: cumulative travel time at each point (`ta[1] = 0`), length n.
- `tra`: total travel time (= `ta[end]`).
"""
function murat_traveling(x, y, z, grid, v, pvel)
    xG, yG, zG = grid.x, grid.y, grid.z
    distlim = 5 * maximum((xG[2] - xG[1], yG[2] - yG[1], zG[2] - zG[1]))
    n = length(x)
    ta = zeros(n)
    ds = sqrt((x[2] - x[1])^2 + (y[2] - y[1])^2 + (z[2] - z[1])^2)
    if ds <= distlim
        # As in MATLAB: the first segment uses the mean slowness, the following
        # ones the inverse of the mean velocity of their end points
        tra = ds * ((1 / v[1] + 1 / v[2]) / 2)
        ta[2] = tra
        for i in 3:n
            ds = sqrt((x[i] - x[i-1])^2 + (y[i] - y[i-1])^2 + (z[i] - z[i-1])^2)
            tra += ds / ((v[i] + v[i-1]) / 2)
            ta[i] = tra
        end
    else
        # km/s -> m/s
        pvelm = pvel .* 1000
        tra = 0.0
        for i in 2:n
            xd, yd, zd = x[i] - x[i-1], y[i] - y[i-1], z[i] - z[i-1]
            ds = sqrt(xd^2 + yd^2 + zd^2)
            # Direction cosines of the segment, split into `is` equal pieces
            xc, yc, zc = xd / ds, yd / ds, zd / ds
            is = trunc(Int, ds / distlim)
            ds /= is
            for j in 1:is
                # Velocity at the midpoint of piece j
                vp = murat_velocity(x[i-1] + (j - 0.5) * ds * xc, y[i-1] + (j - 0.5) * ds * yc,
                                    z[i-1] + (j - 0.5) * ds * zc, grid, pvelm)
                tra += ds / vp
            end
            ta[i] = tra
        end
    end
    return ta, tra
end

"""
    murat_bending!(x, y, z, grid, v, pvel) -> (x, y, z, v)

One pseudo-bending pass (Um and Thurber, 1987; MATLAB header cites Block, 1991):
each interior ray point is moved to the midpoint of its two neighbours and then
displaced along the component of the velocity gradient normal to the ray, by the
distance that minimizes the travel time of the local three-point segment. This
bends the ray towards faster regions. Modifies `x`, `y`, `z`, `v` in place.
MATLAB: `Murat_bending.m` (its convergence-enhancement factor `xfac = 1` has no
effect and is omitted).

# Arguments
- `x`, `y`, `z`: ray point coordinates in m (end points are kept fixed).
- `grid`, `pvel`: propagation grid and velocity model (km/s).
- `v`: velocity at each ray point (km/s), updated for the moved points.

# Returns
- The same (modified) `x`, `y`, `z`, `v`.
"""
function murat_bending!(x, y, z, grid, v, pvel)
    n = length(x)
    for k in 2:n-1
        kk, kkk = k - 1, k + 1
        # Chord between the neighbours: squared length dn, unit tangent (rdx, rdy, rdz)
        dx, dy, dz = x[kkk] - x[kk], y[kkk] - y[kk], z[kkk] - z[kk]
        dn = dx^2 + dy^2 + dz^2
        ddn = sqrt(dn)
        rdx, rdy, rdz = dx / ddn, dy / ddn, dz / ddn
        # Midpoint of the chord, its velocity and velocity gradient
        xk, yk, zk = x[kk] + 0.5dx, y[kk] + 0.5dy, z[kk] + 0.5dz
        vk = murat_velocity(xk, yk, zk, grid, pvel)
        vx, vy, vz = murat_veldicity(xk, yk, zk, grid, pvel)
        # Gradient component perpendicular to the tangent (bending direction)
        vrd = vx * rdx + vy * rdy + vz * rdz
        rvx, rvy, rvz = vx - vrd * rdx, vy - vrd * rdy, vz - vrd * rdz
        rvs = sqrt(rvx^2 + rvy^2 + rvz^2)
        if rvs == 0
            # No normal gradient: straight path through the midpoint
            x[k], y[k], z[k], v[k] = xk, yk, zk, vk
        else
            rvx, rvy, rvz = rvx / rvs, rvy / rvs, rvz / rvs
            # Optimal displacement rtemp along the normal (Um & Thurber formula),
            # c = mean slowness of the two neighbours
            c = (1 / v[kkk] + 1 / v[kk]) * 0.5
            c1 = (c * vk + 1) / (4 * c * rvs)
            rtemp = -c1 + sqrt(c1^2 + dn / (8 * c * vk))
            x[k], y[k], z[k] = xk + rvx * rtemp, yk + rvy * rtemp, zk + rvz * rtemp
            v[k] = murat_velocity(x[k], y[k], z[k], grid, pvel)
        end
    end
    return x, y, z, v
end

"""
    murat_tracing(ray, grid, pvel) -> rma (n × 5)

Minimum-time path between source (`ray[:,1]`) and receiver (`ray[:,2]`) by
pseudo-bending. Starts from the straight segment (3 points), bends it repeatedly
(`murat_bending!`) and, whenever the travel time stops decreasing by more than the
time threshold, doubles the number of points (inserting midpoints), until both the
travel time and the point spacing converge (or 100 iterations).
MATLAB: `Murat_tracing.m`.

# Arguments
- `ray`: 3 × 2 matrix, columns = source and receiver `[x; y; z]` in m.
- `grid`: propagation grid (see `murat_cornering`).
- `pvel`: propagation velocity model (ny, nx, nz), km/s.

# Returns
- `rma`: n × 5 matrix `[index x y z time]` per ray point; x, y, z in m, time =
  cumulative travel time from `murat_traveling` divided by 1000 (s when the short
  branch of `murat_traveling` is used).
"""
function murat_tracing(ray, grid, pvel)
    maxit = 100
    # Maximum number of ray points (30 × 10000 in MATLAB)
    max3 = 3 * 10000
    xG, yG, zG = grid.x, grid.y, grid.z
    # dconv: diagonal of one propagation cell (m); tlim: time to cross it (s) at
    # the mean model velocity (pvel in km/s -> m/s)
    dconv = sqrt((xG[2] - xG[1])^2 + (yG[2] - yG[1])^2 + (zG[2] - zG[1])^2)
    tlim = dconv / mean(pvel .* 1000)
    # Initial ray: source, midpoint, receiver
    x = [ray[1, 1], (ray[1, 1] + ray[1, 2]) / 2, ray[1, 2]]
    y = [ray[2, 1], (ray[2, 1] + ray[2, 2]) / 2, ray[2, 2]]
    z = [ray[3, 1], (ray[3, 1] + ray[3, 2]) / 2, ray[3, 2]]
    n = 3
    v = [murat_velocity(x[i], y[i], z[i], grid, pvel) for i in 1:n]
    ta, tra = murat_traveling(x, y, z, grid, v, pvel)
    ttsave = 0.0
    for _ in 1:maxit
        prevT = tra
        murat_bending!(x, y, z, grid, v, pvel)
        ta, tra = murat_traveling(x, y, z, grid, v, pvel)
        # Bending no longer reduces the travel time significantly
        if prevT - tra < tlim
            # Converged if the time is stable since the last refinement and the
            # point spacing is below a tenth of a cell diagonal
            dtemp = sqrt((x[2] - x[1])^2 + (y[2] - y[1])^2 + (z[2] - z[1])^2)
            if abs(ttsave - tra) < tlim && dtemp < dconv / 10
                break
            end
            # Otherwise double the number of segments
            ttsave = tra
            nnew = 2n - 1
            nnew > max3 && error("Too many points, reduce them or increase max3")
            # Old points at odd indices, midpoints of consecutive points at even ones
            refine(a) = [isodd(i) ? a[(i+1)÷2] : 0.5 * (a[i÷2] + a[i÷2+1]) for i in 1:nnew]
            x, y, z = refine(x), refine(y), refine(z)
            n = nnew
            v = [murat_velocity(x[i], y[i], z[i], grid, pvel) for i in 1:n]
            ta = zeros(n)
        end
    end
    # As MATLAB: travel time divided by 1000 (ms -> s for the short branch)
    return hcat(collect(1.0:n), x, y, z, ta ./ 1000)
end

"""
    murat_segments(modv, rma) -> (blocch, lunpar, luntot, s, rayCrossing)

Length of the ray inside each cell of the inversion grid. Each polyline segment
of the traced ray is split at every crossed x/y/z cell boundary and each piece is
assigned to the cell containing its midpoint. These lengths fill the peak-delay
and Q inversion matrices in `murat_rays`.
MATLAB: `Murat_segments.m`. Note the output order differs from MATLAB
(`[lunpar, blocch, ...]` there).

As in MuRAT, a sub-segment is assigned to the cell following the one whose
node range contains its midpoint (`block = index + 1`).

# Arguments
- `modv`: inversion model, nModel × 4 `[x y z v]` (m, m, m altitude, km/s), in MuRAT
  list order (x slowest, z fastest, z decreasing).
- `rma`: traced ray from `murat_tracing`, n × 5 `[index x y z time]`.

# Returns
- `blocch`: indices (rows of `modv`) of the crossed cells.
- `lunpar`: ray length in each crossed cell (m).
- `luntot`: total ray length inside the model (m), `sum(lunpar)`.
- `s`: slowness of each crossed cell, `1 ./ modv[blocch, 4]` (s/km).
- `rayCrossing`: number of sub-segments in each of the nModel cells (hit count).
"""
function murat_segments(modv, rma)
    nBlocks = size(modv, 1)
    lunByBlock = zeros(nBlocks)
    rayCrossing = zeros(nBlocks)
    xv, yv, zv = modv[:, 1], modv[:, 2], modv[:, 3]
    lgx, lgy, lvg = maximum(xv), maximum(yv), minimum(zv)
    lgx1, lgy1, lvg1 = minimum(xv), minimum(yv), maximum(zv)
    # Grid steps from the list order: z changes every row, y every nz rows
    # (passoy = nz), x every ny*nz rows (passox); dsz is negative
    passox = findfirst(!=(xv[1]), xv) - 1
    passoy = findfirst(!=(yv[1]), yv) - 1
    dsx = modv[1+passox, 1] - modv[1, 1]
    dsy = modv[1+passoy, 2] - modv[1, 2]
    dsz = modv[2, 3] - modv[1, 3]
    # Cell boundary planes along x, y (both extended to the largest grid extent,
    # as MATLAB) and z (from the top down)
    mn1 = maximum((lgx - lgx1, lgy - lgy1, lvg - lvg1))
    BLx = matlab_colon(lgx1, dsx, lgx1 + mn1)
    BLy = matlab_colon(lgy1, dsy, lgy1 + mn1)
    BLv = matlab_colon(lvg1, dsz, lvg)

    # Keep only ray points with x > 0 and y > 0 (as MATLAB), then drop repeated points
    fallr = findall((rma[:, 2] .> 0) .& (rma[:, 3] .> 0))
    if length(fallr) < 2
        return Int[], Float64[], 0.0, Float64[], rayCrossing
    end
    pts = rma[fallr, 2:4]
    step = [0.0; sqrt.(sum(diff(pts, dims = 1) .^ 2, dims = 2))[:]]
    pts = pts[[true; step[2:end] .> 0], :]

    for ir in 1:size(pts, 1)-1
        # Segment p0 -> p0 + dp, parametrized by t in [0, 1]
        p0 = pts[ir, :]
        dp = pts[ir+1, :] .- p0
        segLength = sqrt(sum(abs2, dp))
        segLength == 0 && continue
        # Parameters t where the segment crosses a boundary plane, plus the end points
        tB = [0.0, 1.0]
        for (BL, d) in ((BLx, 1), (BLy, 2), (BLv, 3))
            if dp[d] != 0
                tt = (BL .- p0[d]) ./ dp[d]
                append!(tB, tt[(tt .> 0) .& (tt .< 1)])
            end
        end
        sort!(tB)
        # Merge crossings closer than 1e-10 (e.g. through an edge or corner)
        tB = tB[[true; diff(tB) .> 1e-10]]
        for it in 1:length(tB)-1
            t0, t1 = tB[it], tB[it+1]
            t1 <= t0 && continue
            # First node whose cell range [x, x+dsx) × [y, y+dsy) × [z+dsz, z)
            # contains the piece midpoint
            pm = p0 .+ (t0 + t1) / 2 .* dp
            idx = findfirst(b -> xv[b] <= pm[1] < xv[b] + dsx && yv[b] <= pm[2] < yv[b] + dsy &&
                                 zv[b] + dsz <= pm[3] < zv[b], 1:nBlocks)
            idx === nothing && continue
            # MuRAT quirk kept on purpose: length assigned to the next row of modv
            block = idx + 1
            (block > nBlocks || block < 1) && continue
            lunByBlock[block] += segLength * (t1 - t0)
            rayCrossing[block] += 1
        end
    end
    blocch = findall(>(0), lunByBlock)
    lunpar = lunByBlock[blocch]
    s = isempty(blocch) ? Float64[] : 1 ./ modv[blocch, 4]
    return blocch, lunpar, sum(lunpar), s, rayCrossing
end

"""
    murat_rays(modv, grid, pvel, locM) -> (Apd, AQ, luntot, rplot, rayCrossing)

Traces the ray of one source-station pair (`murat_tracing`) and measures its
length in each inversion cell (`murat_segments`) to build that pair's rows of the
peak-delay and Q (direct-wave attenuation) inversion matrices. Called once per
waveform in `murat_data` (`data.jl`).
MATLAB: `Murat_rays.m`.

# Arguments
- `modv`: inversion model, nModel × 4 `[x y z v]` (m, km/s).
- `grid`: propagation grid (see `murat_cornering`).
- `pvel`: propagation velocity model (ny, nx, nz), km/s.
- `locM`: `[xs ys zs xr yr zr]`, source then receiver coordinates in m.

# Returns
- `Apd`: peak-delay matrix row, ray length per cell in km (length nModel).
- `AQ`: Q matrix row, `-length (km) × slowness (s/km)` per cell = minus the travel
  time spent in each cell (s) (length nModel).
- `luntot`: total ray length inside the model (m).
- `rplot`: 100 × 5 subsample of the ray `[index x y z time]` for plotting.
- `rayCrossing`: hit count per cell (length nModel).
"""
function murat_rays(modv, grid, pvel, locM)
    nCells = size(modv, 1)
    Apd = zeros(nCells)
    AQ = zeros(nCells)
    # Columns: source, receiver
    ray = reshape(locM, 3, 2)
    rma = murat_tracing(ray, grid, pvel)
    blocch, lunpar, luntot, s, rayCrossing = murat_segments(modv, rma)
    # Lengths m -> km
    Apd[blocch] .= lunpar ./ 1000
    AQ[blocch] .= -lunpar ./ 1000 .* s
    # 100 equally spaced ray points (MATLAB round: ties away from zero); points are
    # repeated if the ray has fewer than 100
    lrma = rma[end, 1]
    idx = round.(Int, range(1, max(1, lrma), length = 100), RoundNearestTiesAway)
    return Apd, AQ, luntot, rma[idx, :], rayCrossing
end
