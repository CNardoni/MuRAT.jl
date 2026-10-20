# Figures of a MuRAT run with CairoMakie (loaded with `using CairoMakie`).
# Package extension: Julia loads this module automatically only when both MuRAT
# and CairoMakie are loaded; it adds the method `MuRAT.murat_plot(m)`, the
# translation of `Murat_plot.m`. Figures are PNG files in the same folders as
# MATLAB, but drawn with Makie (same content, not copies of the MATLAB figures).
module MuratMakieExt

using MuRAT
using MuRAT: Murat, outdir, fmt_g, to_xyz_array, murat_lsqlin_qmean, georeckon, fold_field
using CairoMakie
using Statistics
using Printf

"""
    fcname(f) -> String

Frequency as used in file names, e.g. `1.5` -> `"1_5"`, `6.0` -> `"6"` (as
MATLAB `strrep(num2str(f),'.','_')`). No MATLAB counterpart (Julia helper).
"""
fcname(f) = replace(string(fmt_g(f)), "." => "_")


# Approximate km per degree, for the panel aspect ratios only
const KM_PER_DEG_LAT = 111.0
km_per_deg_lon(lat) = 111.32 * cosd(lat)

"""
    section_planes(m::Murat) -> (iz, iy, ix)

Indices of the grid planes nearest to `m.input.sections = [lat, lon, altitude (m)]`:
altitude plane of the horizontal section, latitude plane of the west-east
section, longitude plane of the south-north section. Sections are drawn on these
planes, so the cells keep their own values (no averaging of two planes).
"""
function section_planes(m::Murat)
    g, sec = m.grid, m.input.sections
    near(c, v) = argmin(abs.(c .- v))
    return near(g.z, sec[3]), near(g.y, sec[1]), near(g.x, sec[2])
end

"""
    plane_positions(m::Murat, planes) -> (altitude_m, lat, lon)

Positions of the grid planes `planes = (iz, iy, ix)` (see `section_planes`).
"""
plane_positions(m::Murat, planes) = (m.grid.z[planes[1]], m.grid.y[planes[2]], m.grid.x[planes[3]])

"""
    slice_at(A, coords, value, dim) -> 2D array

Section of the 3D field `A` (nx × ny × nz) at `coords = value` along dimension
`dim`, linearly interpolated between the two nearest grid planes (`value` is
clamped to the grid; `coords` may be decreasing, e.g. altitude).
"""
function slice_at(A, coords, value, dim)
    n = length(coords)
    if coords[end] < coords[1]
        coords, A = reverse(coords), reverse(A, dims = dim)
    end
    v = clamp(value, coords[1], coords[end])
    i = clamp(searchsortedlast(coords, v), 1, n - 1)
    f = (v - coords[i]) / (coords[i+1] - coords[i])
    return (1 - f) .* selectdim(A, dim, i) .+ f .* selectdim(A, dim, i + 1)
end

"""
    refine(c, f) -> Vector

Grid vector `c` with `f` times smaller steps (same first and last values).
"""
refine(c, f) = collect(range(c[1], c[end], length = (length(c) - 1) * f + 1))

"""
    bilinear(a, b, M, aq, bq) -> Matrix

Bilinear interpolation of `M` (length(a) × length(b), on grid vectors `a`, `b`,
increasing or decreasing) at all the points of the vectors `aq` × `bq`.
"""
function bilinear(a, b, M, aq, bq)
    function locate(c, q)
        rev = c[end] < c[1]
        cc = rev ? reverse(c) : c
        i = clamp(searchsortedlast(cc, q), 1, length(cc) - 1)
        f = clamp((q - cc[i]) / (cc[i+1] - cc[i]), 0, 1)
        # Index and weight in the original (possibly decreasing) order
        return rev ? (length(c) - i, 1 - f) : (i, f)
    end
    out = zeros(length(aq), length(bq))
    for (j, bv) in enumerate(bq), (i, av) in enumerate(aq)
        ia, fa = locate(a, av)
        ib, fb = locate(b, bv)
        out[i, j] = (1 - fa) * (1 - fb) * M[ia, ib] + fa * (1 - fb) * M[ia+1, ib] +
                    (1 - fa) * fb * M[ia, ib+1] + fa * fb * M[ia+1, ib+1]
    end
    return out
end

"""
    sections3(m::Murat, v; planes = section_planes(m), interp = 1)
        -> (x, y, d, H, WE, SN, at)

Horizontal, west-east and south-north sections of a field in MuRAT list order.
With `interp = 1` they are the grid planes nearest to `m.input.sections` (cells
keep their own values); with `interp > 1` they are taken at the exact section
positions (linear between planes) and bilinearly interpolated on a grid `interp`
times finer, as the interpolated slices of `Murat_image3D.m`.
No MATLAB counterpart (Julia helper; MATLAB folds with `Murat_fold` and uses
`slice` in `Murat_image3D.m`).

# Arguments
- `v`: field value per grid cell, MuRAT list order (x slowest, z fastest).
- `planes`: plane indices `(iz, iy, ix)` used when `interp = 1`.
- `interp`: refinement factor (input `plotInterpolation`).

# Returns
- `x`, `y`: longitudes and latitudes (°); `d`: depths (km, positive down,
  increasing) of the plotted nodes.
- `H` (nx × ny), `WE` (nx × nd), `SN` (ny × nd): the three sections.
- `at`: section positions `(altitude (m), lat, lon)`.
"""
function sections3(m::Murat, v; planes = section_planes(m), interp = 1)
    g = m.grid
    x, y, z = g.x, g.y, g.z ./ 1000
    A = to_xyz_array(v, length(x), length(y), length(z))
    if interp <= 1
        iz, iy, ix = planes
        H, WE, SN = A[:, :, iz], A[:, iy, :], A[ix, :, :]
        at = plane_positions(m, planes)
    else
        sec = m.input.sections
        H = slice_at(A, z, sec[3] / 1000, 3)       # (nx, ny)
        WE = slice_at(A, y, sec[1], 2)             # (nx, nz)
        SN = slice_at(A, x, sec[2], 1)             # (ny, nz)
        xr, yr, zr = refine(x, interp), refine(y, interp), refine(z, interp)
        H, WE, SN = bilinear(x, y, H, xr, yr), bilinear(x, z, WE, xr, zr), bilinear(y, z, SN, yr, zr)
        x, y, z = xr, yr, zr
        at = (clamp(sec[3], extrema(g.z)...), clamp(sec[1], extrema(g.y)...),
              clamp(sec[2], extrema(g.x)...))
    end
    # Altitude (km, up) to depth (km, down), sorted by increasing depth
    d = -z
    p = sortperm(d)
    return x, y, d[p], H, WE[:, p], SN[:, p], at
end

"""
    topo_grid(topo) -> (lon, lat, elev) or nothing

Regular topography grid from the `topo` keyword of `murat_plot`: longitudes and
latitudes (°, increasing) and elevation (km, positive up, lon × lat).
Accepts the `GeoData` of GeophysicalModelGenerator `import_topo` (fields `lon`,
`lat`, `depth` in km; also when read back from a JLD2 file without
GeophysicalModelGenerator loaded) or a `NamedTuple` `(lon, lat, elevation)` with
vectors and an elevation matrix in km.
No MATLAB counterpart (Julia addition).
"""
function topo_grid(topo)
    topo === nothing && return nothing
    # GeoUnit values (GeophysicalModelGenerator) keep the numbers in `.val`
    val(a) = hasproperty(a, :val) ? getproperty(a, :val) : a
    if hasproperty(topo, :elevation)
        lon, lat, el = collect(topo.lon), collect(topo.lat), Matrix(topo.elevation)
    else
        LON, LAT, DEP = val(topo.lon), val(topo.lat), val(topo.depth)
        lon, lat, el = LON[:, 1, 1], LAT[1, :, 1], DEP[:, :, 1]
    end
    # Increasing coordinates for the interpolation
    if lon[end] < lon[1]
        lon, el = reverse(lon), reverse(el, dims = 1)
    end
    if lat[end] < lat[1]
        lat, el = reverse(lat), reverse(el, dims = 2)
    end
    return Float64.(lon), Float64.(lat), Float64.(el)
end

"""
    topo_profile(tg, lons, lats) -> Vector

Elevation (km) at the points `(lons[i], lats[i])`, bilinear interpolation of the
topography grid `tg = (lon, lat, elev)` from `topo_grid` (clamped at its edges).
"""
function topo_profile(tg, lons, lats)
    lon, lat, el = tg
    out = similar(lons, Float64)
    for n in eachindex(lons, lats)
        i = clamp(searchsortedlast(lon, lons[n]), 1, length(lon) - 1)
        j = clamp(searchsortedlast(lat, lats[n]), 1, length(lat) - 1)
        fx = clamp((lons[n] - lon[i]) / (lon[i+1] - lon[i]), 0, 1)
        fy = clamp((lats[n] - lat[j]) / (lat[j+1] - lat[j]), 0, 1)
        out[n] = (1 - fx) * (1 - fy) * el[i, j] + fx * (1 - fy) * el[i+1, j] +
                 (1 - fx) * fy * el[i, j+1] + fx * fy * el[i+1, j+1]
    end
    return out
end

"""
    add_topography!(ax, h, elev, dlims)

Draws a topography profile on a vertical section: elevation `elev` (km) along
the horizontal coordinate `h`, as a black line at depth `-elev` with the air
above it shaded. Returns the depth limits extended to show the highest point.
"""
function add_topography!(ax, h, elev, dlims)
    dtop = min(dlims[1], -maximum(elev) - 0.2)
    band!(ax, h, fill(dtop, length(h)), -elev, color = (:lightblue, 0.35))
    lines!(ax, h, -elev, color = :black, linewidth = 1.5)
    return (dtop, dlims[2])
end

"""
    section_axes(gl, m, title; ve = 2) -> (aH, aWE, aSN, lay)

Layout shared by all maps and the ray figures, in the grid layout `gl`: the
horizontal panel on the left (lon-lat with the aspect of the real distances), the
west-east and south-north sections stacked on the right, depth increasing down,
with vertical exaggeration `ve`.
No MATLAB counterpart (Julia helper; MATLAB draws 3D slices).

# Arguments
- `gl`: a figure position or `GridLayout`; `m`: the run (grid and sections).
- `title`: title of the horizontal panel.
- `ve`: vertical exaggeration of the sections (1 = true scale); raised automatically
  when the sections would be more than 4 times wider than tall.
  Panel sizes follow the grid shape (the figure is then trimmed with
  `resize_to_layout!`).
- `at`: section positions `(altitude (m), lat, lon)` shown in the titles.

# Returns
- the three Makie axes and their `GridLayout` (for the colour bar).
"""
function section_axes(gl, m, title; ve = 2, at = plane_positions(m, section_planes(m)))
    g = m.grid
    x, y = g.x, g.y
    d = sort(-g.z ./ 1000)
    latm = (y[1] + y[end]) / 2
    wx = abs(x[end] - x[1]) * km_per_deg_lon(latm)       # km west-east
    wy = abs(y[end] - y[1]) * KM_PER_DEG_LAT             # km south-north
    wz = (d[end] - d[1]) + 1.5                           # km depth (+ topography)
    lay = GridLayout(gl)
    aH = Axis(lay[1:2, 1], title = title, xlabel = "Longitude (°)", ylabel = "Latitude (°)",
              aspect = AxisAspect(wx / wy))
    aWE = Axis(lay[1, 2], title = @sprintf("West-East (lat %.3f°)", at[2]),
               xlabel = "Longitude (°)", ylabel = "Depth (km)", yreversed = true)
    aSN = Axis(lay[2, 2], title = @sprintf("South-North (lon %.3f°)", at[3]),
               xlabel = "Latitude (°)", ylabel = "Depth (km)", yreversed = true)
    # Panel sizes (px) from the grid shape: the two sections share one width and
    # together fill the height of the map, so no space is left empty whatever the
    # grid (wide and shallow, or small and deep)
    H, gap = 560.0, 75.0                     # section heights; labels between them
    # Vertical exaggeration: `ve`, raised for wide and shallow grids so that the
    # sections are at most 4 times wider than tall
    ve = max(ve, max(wx, wy) / (4 * wz))
    rWE, rSN = wx / (wz * ve), wy / (wz * ve)            # width/height of the sections
    W = H / (1 / rWE + 1 / rSN)
    rowsize!(lay, 1, Fixed(W / rWE))
    rowsize!(lay, 2, Fixed(W / rSN))
    colsize!(lay, 1, Fixed((H + gap) * wx / wy))         # map with real proportions
    colsize!(lay, 2, Fixed(W))
    return aH, aWE, aSN, lay
end

"""
    finish_sections!(aH, aWE, aSN, m, topo)

Common limits of the three panels (grid extent; depth from the grid top, or the
highest topography, to the grid bottom) and, if `topo` is given, the topography
profiles along the two vertical sections.
"""
function finish_sections!(aH, aWE, aSN, m, topo; at = plane_positions(m, section_planes(m)))
    g = m.grid
    xl, yl = extrema(g.x), extrema(g.y)
    dl = extrema(-g.z ./ 1000)
    dWE = dSN = dl
    tg = topo_grid(topo)
    if tg !== nothing
        hx = range(xl..., length = 400)
        hy = range(yl..., length = 400)
        dWE = add_topography!(aWE, hx, topo_profile(tg, hx, fill(at[2], 400)), dl)
        dSN = add_topography!(aSN, hy, topo_profile(tg, fill(at[3], 400), hy), dl)
    end
    limits!(aH, xl..., yl...)
    # Depth limits given bottom first: limits! would undo `yreversed`
    xlims!(aWE, xl...)
    ylims!(aWE, dWE[2], dWE[1])
    xlims!(aSN, yl...)
    ylims!(aSN, dSN[2], dSN[1])
    return nothing
end

"""
    draw_field!(gl, m, v, title; colormap, colorrange = nothing, events = nothing,
                topo = nothing, label = "") -> heatmap

One field on the three panels of `section_axes`, with a colour bar: flat coloured
cells centred on the grid nodes (`interp = 1`), or on a grid `interp` times finer
after bilinear interpolation (`interp > 1`, see `sections3`).

# Arguments
- `v`: field per grid cell in MuRAT list order (e.g. `m.Qc.modvQc[:, 4, k]`).
- `colorrange`: `(min, max)`, or `nothing` for the range of the shown values.
- `events`: `[evla evlo evdepth stla stlo stdepth]` (km) to plot on the horizontal
  panel (events as circles, stations as triangles), or `nothing`.
- `topo`: topography for the vertical sections (see `topo_grid`), or `nothing`.
- `label`: colour-bar label.
- `planes`: grid-plane indices `(iz, iy, ix)` of the sections when `interp = 1`
  (default: nearest to `m.input.sections`, see `section_planes`).
- `interp`: interpolation factor (1 = flat cells on the grid).
"""
function draw_field!(gl, m, v, title; colormap = :inferno, colorrange = nothing,
                     events = nothing, topo = nothing, label = "",
                     planes = section_planes(m), interp = 1)
    x, y, d, H, WE, SN, at = sections3(m, v; planes = planes, interp = interp)
    cr = colorrange === nothing ?
         extrema(filter(isfinite, vcat(vec(H), vec(WE), vec(SN)))) : colorrange
    # Avoid an empty colour range for a constant field
    cr[1] == cr[2] && (cr = (cr[1] - 1e-12, cr[2] + 1e-12))
    aH, aWE, aSN, lay = section_axes(gl, m, @sprintf("%s – depth %.1f km", title, -at[1] / 1000);
                                     at = at)
    hm = heatmap!(aH, x, y, H, colormap = colormap, colorrange = cr)
    heatmap!(aWE, x, d, WE, colormap = colormap, colorrange = cr)
    heatmap!(aSN, y, d, SN, colormap = colormap, colorrange = cr)
    # Events and stations on the horizontal panel
    if events !== nothing
        scatter!(aH, events[:, 2], events[:, 1], color = :cyan, strokecolor = :blue,
                 strokewidth = 1, markersize = 6)
        scatter!(aH, events[:, 5], events[:, 4], marker = :utriangle, color = :gray,
                 strokecolor = :magenta, strokewidth = 1, markersize = 10)
    end
    finish_sections!(aH, aWE, aSN, m, topo; at = at)
    Colorbar(lay[1:2, 3], hm, label = label)
    return hm
end

"""
    map_figure(m::Murat, v, title; colormap = :inferno, colorrange = nothing,
               events = nothing, topo = nothing, label = "") -> Figure

Map of a 3D field: horizontal, west-east and south-north sections as flat cells,
depth positive down, topography on the vertical sections when `topo` is given.
Used for the result maps (peak delay, Qc, Q), the peak-delay hit map and the
kernels.
MATLAB: `Murat_image3D.m` (also `Murat_hitmap.m` and `Murat_imageKernels.m`);
MATLAB draws 3D slices (interpolated and smoothed for the results), here three
2D panels without interpolation.
"""
function map_figure(m::Murat, v, title; kw...)
    fig = Figure(size = (1600, 640))
    draw_field!(fig[1, 1], m, v, title; kw...)
    resize_to_layout!(fig)
    return fig
end

"""
    test_figure(m::Murat, vin, vout, title; colormap = :bone, topo = nothing, interp = 1) -> Figure

Input (top) and output (bottom) of a synthetic test (checkerboard or spike), each
on the three sections of `draw_field!` with its own colour scale.
MATLAB: the checkerboard and spike figures of `Murat_plot.m`, drawn with
`Murat_image3D_2panels.m`.

# Arguments
- `vin`, `vout`: input and recovered model per cell, MuRAT list order (e.g.
  columns 6 and 7 of `m.Qc.modvQc[:, :, k]`).
- `title`: test name used in the panel titles.
"""
function test_figure(m::Murat, vin, vout, title; colormap = :bone, topo = nothing, interp = 1)
    fig = Figure(size = (1600, 1250))
    draw_field!(fig[1, 1], m, vin, "Input $title"; colormap = colormap, topo = topo, interp = interp)
    draw_field!(fig[2, 1], m, vout, "Output $title"; colormap = colormap, topo = topo, interp = interp)
    resize_to_layout!(fig)
    return fig
end

"""
    rays_figure(m::Murat, mask, title; topo = nothing) -> Figure

Ray paths of the selected traces in map view (lon-lat) and projected on the
west-east (lon-depth) and south-north (lat-depth) planes, with events and
stations on the map and the topography profiles of the sections when `topo` is
given.
MATLAB: `Murat_imageRays.m` (rays converted from meters to lon/lat with
`reckon` on WGS84 and azimuth `atan(x/y)`, kept for compatibility). The
coastline and location panels of MATLAB are not drawn.

# Arguments
- `mask`: Bool vector, traces whose rays are drawn (e.g. `m.Q.retainQ[:, k]`).
- `title`: title of the map panel.
"""
function rays_figure(m::Murat, mask, title; topo = nothing)
    s = m.input
    rp = m.rays.raysPlot
    fig = Figure(size = (1600, 640))
    aH, aWE, aSN, _ = section_axes(fig[1, 1], m, title)
    # All rays in one NaN-separated line per panel: one plot object instead of one per
    # ray (tens of thousands of lines! calls take hours in Makie)
    idx = findall(mask)
    np = size(rp, 1) + 1
    lat, lon, dd = (fill(NaN, np * length(idx)) for _ in 1:3)
    for (j, i) in enumerate(idx)
        o = (j - 1) * np
        # raysPlot[:, :, i] columns: point index, x, y (m from the origin),
        # altitude (m), time; atan(x/y) is NaN at the origin -> azimuth 0
        for n in 1:np-1
            xn, yn = rp[n, 2, i], rp[n, 3, i]
            az = rad2deg(atan(xn / yn))
            lat[o+n], lon[o+n] = georeckon(s.origin[1], s.origin[2], sqrt(xn^2 + yn^2), isnan(az) ? 0.0 : az)
            dd[o+n] = -rp[n, 4, i] / 1000             # depth (km)
        end
    end
    lines!(aH, lon, lat, color = (:black, 0.4), linewidth = 0.6)
    lines!(aWE, lon, dd, color = (:black, 0.4), linewidth = 0.6)
    lines!(aSN, lat, dd, color = (:black, 0.4), linewidth = 0.6)
    # Events (circles) and stations (triangles) from locationsDeg
    # [evla evlo evdp stla stlo stel]
    ev = m.rays.locationsDeg[mask, :]
    scatter!(aH, ev[:, 2], ev[:, 1], color = :cyan, strokecolor = :blue, strokewidth = 1, markersize = 6)
    scatter!(aH, ev[:, 5], ev[:, 4], marker = :utriangle, color = :gray, strokecolor = :magenta,
             strokewidth = 1, markersize = 10)
    finish_sections!(aH, aWE, aSN, m, topo)
    resize_to_layout!(fig)
    return fig
end

"""
    log_ticks(lo, hi) -> (values, labels)

Ticks for a log axis written as plain numbers (e.g. 2, 5, 10, 20 s) instead of
powers of ten: 1-2-5 per decade, with 3 added when there are few ticks.
"""
function log_ticks(lo, hi)
    decades = floor(Int, log10(lo)):ceil(Int, log10(hi))
    pick(mants) = filter(t -> lo <= t <= hi, [m * 10.0^e for e in decades for m in mants])
    v = pick((1, 2, 5))
    length(v) < 4 && (v = pick((1, 2, 3, 5)))
    return v, [@sprintf("%g", t) for t in v]
end

"""
    movmean(y, k) -> Vector

Moving mean over `k` points, as MATLAB `movmean(y, k)`: centred window (for even
`k`, k/2 points before and k/2 - 1 after), shrinking at the ends.
"""
function movmean(y, k)
    n = length(y)
    c = [0.0; cumsum(y)]
    b, a = k ÷ 2, iseven(k) ? k ÷ 2 - 1 : k ÷ 2
    return [(c[min(n, i + a)+1] - c[max(1, i - b)]) / (min(n, i + a) - max(1, i - b) + 1)
            for i in 1:n]
end

"""
    MuRAT.murat_plot(m::Murat; topo = nothing) -> Murat

Figures of a run, saved as PNG under `<label>/` (`Rays`, `Tests`, `Results`,
`Checkerboard`, `Spike`, `Kernels`): for each frequency the rays, the data
checks (Qc vs ray length, peak delay vs travel time, coda normalization,
L-curves), the peak-delay, Qc and Q maps, the checkerboard and spike tests;
then the kernel of the first trace and Qc⁻¹ against frequency.
Last stage of a run, after `murat_inversion` and `murat_summary` (which gives
`m.Qc.averageQcFrequency` and `m.Qc.frequencyFit`). Each group is drawn when
its input flag (`PlotRays`, `PlotTests`, `PlotResults`, `PlotCheckers`,
`PlotSpikes`) is 1.
MATLAB: `Murat_plot.m`. Differences: MATLAB draws rays and spikes when
`PlotRays == 0` / `PlotSpikes == 0`; the parameter-separation map (disabled in
MATLAB) and the `.mat`/`.fig` files are not produced; the Qc-frequency average
is computed in `murat_summary` instead of here; maps are flat cells on three 2D
panels (depth positive down) instead of interpolated 3D slices.

# Arguments
- `m`: run after `murat_inversion` (and `murat_summary` for the Qc-frequency
  figure).
- `topo`: optional topography drawn on all vertical sections (maps, tests, rays):
  the `GeoData` of GeophysicalModelGenerator `import_topo`, e.g.
  `import_topo(lon = [-116.6, -115.1], lat = [32.4, 33.6], file = "@earth_relief_15s")`
  (needs `using GeophysicalModelGenerator`), or `(lon, lat, elevation)` in
  degrees and km (see `topo_grid`).

# Returns
- `m`, unchanged.
"""
function MuRAT.murat_plot(m::Murat; topo = nothing)
    s, g = m.input, m.grid
    cf = s.centralFrequency
    evst = m.rays.locationsDeg
    # Event/station table with altitudes in km: [evla evlo -evdp stla stlo stel]
    evst_km = hcat(evst[:, 1:2], -evst[:, 3] ./ 1000, evst[:, 4:5], evst[:, 6] ./ 1000)
    # Save `fig` as <label>/<parts...>.png
    saveto(fig, parts...) = (save(outdir(m, parts...) * ".png", fig); nothing)

    # Interpolation factor of the result and checkerboard maps (1 = flat cells)
    ip = s.plotInterpolation
    for (k, f) in enumerate(cf)
        fc = fcname(f)
        rtpd, rtQc, rtQ = m.PD.retainPeakDelay[:, k], m.Qc.retainQc[:, k], m.Q.retainQ[:, k]

        if s.PlotRays == 1
            saveto(rays_figure(m, rtpd, "Rays peak delay $f Hz"; topo = topo), "Rays",
                   "Rays_PeakDelay_$(fc)_Hz")
            saveto(rays_figure(m, rtQ, "Rays Q $f Hz"; topo = topo), "Rays", "Rays_Q_$(fc)_Hz")
        end

        if s.PlotTests == 1
            # Qc against ray length, as Murat_imageCheckQc.m: only the points within
            # mean ± 1 std, sorted by ray length, with their 500-point moving mean
            Qm = m.Qc.inverseQc[rtQc, k]
            l = m.rays.totalLengthRay[rtQc] ./ 1000
            mu, sd = mean(Qm), std(Qm)
            p = sortperm(l)
            l, Qm = l[p], Qm[p]
            inside = (Qm .> mu - sd) .& (Qm .< mu + sd)
            fig = Figure(size = (900, 650))
            ax = Axis(fig[1, 1], title = "Dependence of Qc⁻¹ on ray lengths ($f Hz) – misfit reduction " *
                                         @sprintf("%.3g", m.Qc.residualRedQc[k]),
                      xlabel = "Ray length (km)", ylabel = "Qc⁻¹")
            scatter!(ax, l[inside], Qm[inside], color = (:gray, 0.6), markersize = 5,
                     label = "Qc⁻¹ within mean ± 1 std")
            hlines!(ax, [mu], color = :red, linewidth = 2, label = @sprintf("<Qc> = %.1f", 1 / mu))
            hlines!(ax, [mu - sd, mu + sd], color = :red, linestyle = :dash, linewidth = 2,
                    label = "mean ± 1 std")
            any(inside) && lines!(ax, l[inside], movmean(Qm[inside], 500), color = :darkred,
                                  linewidth = 3, label = "Moving mean (500 points)")
            Legend(fig[2, 1], ax, orientation = :horizontal, framevisible = false)
            saveto(fig, "Tests", "Qc", "Qc_analysis_$(fc)_Hz")

            # Peak delay against travel time, with the robust fit
            # log10(PD) = fr[1] log10(t) + fr[2] from the data selection
            t0 = m.rays.travelTime[rtpd]
            pd = m.PD.peakDelay[rtpd, k]
            fr = m.PD.fitrobust[:, k]
            # Log axes labelled with the values in seconds (as Murat_imageCheckPeakDelay.m)
            fig = Figure(size = (900, 600))
            ax = Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = "Travel time (s)", ylabel = "Peak delay (s)",
                      title = @sprintf("Peak delays vs travel times: A(f) = %.3g, B(f) = %.3g", fr[2], fr[1]),
                      xticks = log_ticks(extrema(t0)...), yticks = log_ticks(extrema(pd)...))
            scatter!(ax, t0, pd, color = :black, markersize = 6)
            ts = sort(t0)
            lines!(ax, ts, 10 .^ (fr[1] .* log10.(ts) .+ fr[2]), color = :red, linewidth = 2)
            saveto(fig, "Tests", "PeakDelay", "PD_analysis_$(fc)_Hz")

            # Coda normalization: corrected log-energy ratios (d0) against
            # travel time with the average Q⁻¹ line, and direct energy
            # (direct/coda over coda/noise ratio) against travel time.
            # The ± in the title is Q0[2, 2], a variance, as in MATLAB
            mask = rtQ .& rtQc
            if sum(mask) > 2
                t0q = m.rays.travelTime[mask]
                tCm = m.Qc.tCoda[mask, k]
                d0, Q0, _ = murat_lsqlin_qmean(f, m.rays.totalLengthRay[mask] ./ 1000, t0q,
                                               m.Q.energyRatioBodyCoda[mask, k], m.Qc.inverseQc[mask, k],
                                               tCm, tCm .+ s.codaWindow)
                Ed = m.Q.energyRatioBodyCoda[mask, k] ./ m.Q.energyRatioCodaNoise[mask, k]
                fig = Figure(size = (900, 800))
                ax = Axis(fig[1, 1], xlabel = "Travel time (s)", ylabel = "Corrected log-energy ratio",
                          title = @sprintf("Coda normalization (%g Hz): <Q⁻¹> = %.3g ± %.2g", f, Q0[2, 1], Q0[2, 2]))
                scatter!(ax, t0q, d0, color = :black, markersize = 6)
                ts = sort(t0q)
                lines!(ax, ts, -ts .* Q0[2, 1], color = :red, linewidth = 2)
                ax2 = Axis(fig[2, 1], xscale = log10, xlabel = "Travel time (s)", ylabel = "Log direct energy",
                           title = "Log. direct energy vs travel time")
                scatter!(ax2, t0q, log.(Ed), color = :black, markersize = 6)
                saveto(fig, "Tests", "Q", "CN_analysis_$(fc)_Hz")
            end

            # L-curves
            for name in ("Qc", "Q")
                key = "$(name)_$(fc)_Hz"
                haskey(m.Qc.lcurves, key) || continue
                rho, eta, reg = m.Qc.lcurves[key]
                fig = Figure(size = (800, 600))
                ax = Axis(fig[1, 1], xscale = log10, yscale = log10, title = "L-curve $name $f Hz",
                          xlabel = "Residual norm ||A x - b||", ylabel = "Solution norm ||x||")
                lines!(ax, rho, eta, color = :black)
                # Label about 20 points with their regularization parameter
                ni = max(1, length(rho) ÷ 20)
                idx = ni:ni:length(rho)
                scatter!(ax, rho[idx], eta[idx], marker = :x, color = :red)
                text!(ax, rho[idx], eta[idx], text = [@sprintf("%.2g", r) for r in reg[idx]], fontsize = 9)
                saveto(fig, "Tests", "LCurve", "L-curve_$(name)_$(fc)_Hz")
            end
        end

        if s.PlotResults == 1
            # Peak-delay map: as MATLAB, cells with |ΔlogPD| below 5% of the
            # maximum positive (negative) value are set to zero; symmetric scale
            mpd = m.PD.modvPeakDelay[:, 4, k]
            pos, neg = filter(>=(0), mpd), filter(<(0), mpd)
            thrP = isempty(pos) ? 0.0 : maximum(pos) / 100 * 5
            thrN = isempty(neg) ? 0.0 : minimum(neg) / 100 * 5
            mpdRed = mpd .* ((mpd .> thrP) .| (mpd .< thrN))
            lim = maximum(abs.(mpdRed))
            saveto(map_figure(m, mpdRed, "Peak-delay variations $f Hz"; colormap = Reverse(:RdBu),
                              colorrange = (-lim - 1e-12, lim + 1e-12), events = evst_km[rtpd, :],
                              topo = topo, label = "Δ log₁₀ peak delay", interp = ip),
                   "Results", "PeakDelay", "Peak-Delay-3D_$(fc)_Hz")
            saveto(map_figure(m, m.Qc.modvQc[:, 4, k], "Coda attenuation $f Hz"; colormap = :hot,
                              events = evst_km[rtQc, :], topo = topo, label = "Qc⁻¹",
                              interp = ip),
                   "Results", "Qc", "Qc-3D_$(fc)_Hz")
            saveto(map_figure(m, m.Q.modvQ[:, 4, k], "Total attenuation $f Hz"; colormap = :hot,
                              events = evst_km[rtQ, :], topo = topo, label = "Q⁻¹ variations",
                              interp = ip),
                   "Results", "Q", "Q-3D_$(fc)_Hz")
        end

        # Hit map (col 5 of modvPeakDelay) and checkerboard input/output
        # (cols 6/7 of modvQc and modvQ)
        if s.PlotCheckers == 1
            # Hit map saturated at 50 rays, hot colour map (as Murat_hitmap.m)
            hits = m.PD.modvPeakDelay[:, 5, k]
            saveto(map_figure(m, hits, "Peak-delay hit count $f Hz"; colormap = :hot,
                              colorrange = (minimum(hits), 50), topo = topo, label = "Ray count"),
                   "Checkerboard", "PD", "PD-HitMap_$(fc)_Hz")
            saveto(test_figure(m, m.Qc.modvQc[:, 6, k], m.Qc.modvQc[:, 7, k], "checkerboard Qc $f Hz";
                               topo = topo, interp = ip),
                   "Checkerboard", "Qc", "Qc-Checkerboard_$(fc)_Hz")
            saveto(test_figure(m, m.Q.modvQ[:, 6, k], m.Q.modvQ[:, 7, k], "checkerboard Q $f Hz";
                               topo = topo, interp = ip),
                   "Checkerboard", "Q", "Q-Checkerboard_$(fc)_Hz")
        end

        # Spike input/output (cols 8/9), only if a spike was defined
        if s.PlotSpikes == 1 && s.spikeLocationOrigin !== nothing
            saveto(test_figure(m, m.Qc.modvQc[:, 8, k], m.Qc.modvQc[:, 9, k], "spike Qc $f Hz";
                               colormap = :winter, topo = topo),
                   "Spike", "Qc", "Qc-Spike_$(fc)_Hz")
            saveto(test_figure(m, m.Q.modvQ[:, 8, k], m.Q.modvQ[:, 9, k], "spike Q $f Hz";
                               colormap = :hot, topo = topo),
                   "Spike", "Q", "Q-Spike_$(fc)_Hz")
        end
    end

    # Kernels of the first trace on the inversion grid
    fk = m.Qc.firstKernel
    if fk !== nothing
        for (k, f) in enumerate(cf)
            # Qc matrix row of the first trace: kernel on the inversion grid, sum 1
            a = fk.A[:, k]
            kname = fk.type == "differential" ? "Differential kernel" : "Normalized kernel"
            # Sections through the kernel maximum (the input sections may miss it)
            A = to_xyz_array(a, length(g.x), length(g.y), length(g.z))
            ci = argmax(A)
            kplanes = (ci[3], ci[2], ci[1])
            saveto(map_figure(m, a, "$kname, first trace, $f Hz"; colormap = :inferno,
                              topo = topo, label = "Normalized sensitivity", planes = kplanes),
                   "Kernels", "Kernel_$(fcname(f))_Hz")
        end
    end

    # Qc against frequency: weighted mean Qc⁻¹ ± std (av[1, :], av[2, :]) and
    # the power-law fit Qc⁻¹ = Q0 f^n, both from `murat_summary`
    if haskey(m.Qc, :averageQcFrequency) && length(cf) >= 2
        av = m.Qc.averageQcFrequency
        ff = m.Qc.frequencyFit
        fig = Figure(size = (800, 600))
        ax = Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = "Frequency (Hz)", ylabel = "Qc⁻¹",
                  title = @sprintf("Qc⁻¹ = (%.3g ± %.2g) f^(%.3g ± %.2g)", ff.Q0, ff.Q0_uncertainty,
                                   ff.exponent, ff.exponent_uncertainty))
        # Lower error bar capped below the value (log axis cannot show ≤ 0)
        errorbars!(ax, cf, av[1, :], min.(av[2, :], 0.999 .* av[1, :]), av[2, :], whiskerwidth = 8)
        scatter!(ax, cf, av[1, :], color = :black, markersize = 12)
        fs = range(minimum(cf), maximum(cf), length = 50)
        lines!(ax, fs, ff.Q0 .* fs .^ ff.exponent, color = :black, linewidth = 2)
        saveto(fig, "Results", "Qc_vs_frequency")
    end
    return m
end

end # module
