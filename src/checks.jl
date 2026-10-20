# Murat_checks and the velocity models (Murat_modv1D, Murat_modv3D).

# Subfolders of `<label>` created by `murat_checks` (MATLAB `subdirs` in Murat_checks.m,
# plus `temp`, which MATLAB keeps in the working directory)
const OUTPUT_FOLDERS = ["", "Rays", "Kernels", "Tests", "Tests/PeakDelay", "Tests/Qc",
    "Tests/Q", "Tests/LCurve", "Tests/Clustering", "Results", "Results/PeakDelay",
    "Results/Qc", "Results/Q", "Checkerboard", "Checkerboard/Qc",
    "Checkerboard/PD", "Checkerboard/Q", "Spike", "Spike/Qc", "Spike/Q", "TXT", "temp"]

"""
    murat_checks(input) -> Murat

`Murat_checks`: prepares the output folders, reads and checks the SAC
headers, builds the inversion grid and the velocity models.

An existing output folder `label` is moved to `label_old` (replacing a
previous `label_old`).
MATLAB: `Murat_checks.m`. Differences: derived settings are stored in `m.grid` (MATLAB
adds them to `Murat.input`); header times are stored as field names (`"a"`), not as
`"SAChdr.times.a"`; the header table is a tab-separated `Tests/DataHeaders.txt`
(MATLAB: `DataHeaders.xlsx`); a warning is given also when both origin and S times
are missing (MATLAB warns only for flag 1 or 2); LITHO1.0 models are not supported.

# Arguments
- `input`: user settings (`MuratInput`).

# Returns
- a new `Murat` whose `grid` section contains: `originTime`, `PTime`, `STime`,
  `listSac` (SAC paths), `sacHeader` (`SacHeader` per file), `x`, `y` (grid lon/lat
  vectors, degrees), `z` (altitude, m), `gridStepX`, `gridStepY` (m),
  `gridPropagation` (ray-tracing grid vectors, m), `modv`, `modvp` (velocity models as
  n×4 lists `[x y z v]`, m and km/s), `modvPlot` (`modv` with lon/lat), `pvel`
  (propagation model as meshgrid array) and `DDcoordinates` (nodes in degrees).
"""
function murat_checks(input::MuratInput)
    m = Murat(input)
    s = input
    g = m.grid

    # Output folders
    if isdir(s.label)
        old = s.label * "_old"
        isdir(old) && rm(old, recursive = true)
        mv(s.label, old)
        @info "Folder \"$(s.label)\" found: moved to \"$old\""
    end
    for f in OUTPUT_FOLDERS
        mkpath(joinpath(s.label, f))
    end

    # Header times
    g.originTime = s.originTime === nothing ? "o" : s.originTime
    g.PTime = s.PTime
    g.STime = s.STime === nothing ? "t0" : s.STime

    # Data
    g.listSac = list_sac(s.dataDirectory)
    isempty(g.listSac) && error("No SAC files found in $(s.dataDirectory)")
    headers, table, flag = test_data(g.listSac, g.originTime, g.PTime, g.STime)
    g.sacHeader = headers
    flag in (1, 3) && @warn "Missing origin times."
    flag in (2, 3) && @warn "Missing S-wave times."
    write_header_table(joinpath(s.label, "Tests", "DataHeaders.txt"), table)

    # Inversion grid in meters, from the origin: x east, y north (WGS84 distance and
    # azimuth from origin to ending), z altitude (m)
    origin, ending = s.origin, s.ending
    d, az = geodistance(origin[1], origin[2], ending[1], ending[2])
    dist_x = d * sin(deg2rad(az))
    dist_y = d * cos(deg2rad(az))
    xM = collect(range(0, dist_x, length = s.gridLong))
    yM = collect(range(0, dist_y, length = s.gridLat))
    zM = collect(range(origin[3], ending[3], length = s.gridZ))
    modvXYZ = unfold_xyz(xM, yM, zM)

    # Same grid in degrees (x = longitude, y = latitude), used for figures and outputs
    g.x = collect(range(origin[2], ending[2], length = s.gridLong))
    g.y = collect(range(origin[1], ending[1], length = s.gridLat))
    g.z = collect(range(origin[3], ending[3], length = s.gridZ))
    g.gridStepX = xM[2] - xM[1]
    g.gridStepY = yM[2] - yM[1]

    # Velocity models: 1D table (depth km, Vp, Vs) or 3D table (lat, lon, altitude m, v).
    # In 1D the propagation (ray-tracing) grid is the inversion grid; in 3D it has
    # half the inversion grid step.
    velocityModel = joinpath(s.velocityFolder, s.namev)
    if s.availableVelocity == 0
        if s.namev == "LITHO1.0.nc"
            # LITHO1.0 profile at the grid centre; the reader needs NCDatasets
            hasmethod(read_litho1, Tuple{String,Float64,Float64}) ||
                error("LITHO1.0 models need `using NCDatasets` before the run.")
            modvOriginal = read_litho1(velocityModel, mean([s.origin[1], s.ending[1]]),
                                       mean([s.origin[2], s.ending[2]]))
        else
            modvOriginal = read_model_table(velocityModel)
        end
        modv, pvel, modvPlot = murat_modv1d(modvXYZ, modvOriginal, s.POrS, origin)
        g.gridPropagation = (x = xM, y = yM, z = zM)
        g.modv = modv
        g.modvp = modv
        g.modvPlot = modvPlot
        g.pvel = pvel
    elseif s.availableVelocity == 1
        modvOriginal = read_model_table(velocityModel)
        modvP, pvel, modvEqS, modvPlot = murat_modv3d(modvOriginal, origin, ending,
                                                     s.gridLat, s.gridLong, s.gridZ)
        g.gridPropagation = (x = sort(unique(modvP[:, 1])), y = sort(unique(modvP[:, 2])),
                             z = sort(unique(modvP[:, 3]), rev = true))
        g.modv = modvEqS
        g.modvp = modvP
        g.modvPlot = modvPlot
        g.pvel = pvel
    else
        error("availableVelocity must be 0 (1D model) or 1 (3D model)")
    end
    g.DDcoordinates = dd_coordinates(origin, ending, s.gridLat, s.gridLong, s.gridZ)
    return m
end

"""
    write_header_table(path, table)

Writes the table of header values from `test_data` as tab-separated text (one row
per file, missing values left empty).
MATLAB: `writetable(header, 'DataHeaders.xlsx')` in `Murat_checks.m` (Excel file there).

# Arguments
- `path`: output file name.
- `table`: vector of NamedTuples with the same keys (column names).
"""
function write_header_table(path, table)
    open(path, "w") do io
        names = keys(table[1])
        println(io, join(string.(names), '\t'))
        for r in table
            println(io, join([v === missing ? "" : string(v) for v in values(r)], '\t'))
        end
    end
end

"""
    murat_modv1d(modvXYZ, modvOriginal, PorS, origin) -> (modv, pvel, modvPlot)

`Murat_modv1D`: assigns to each grid node the velocity of the closest depth
of a 1D model (columns: depth in km, Vp, Vs).
The same model is used for ray tracing and for the inversion.
MATLAB: `Murat_modv1D.m`.

# Arguments
- `modvXYZ`: grid nodes, n×3 `[x y z]` in m (MuRAT list order, x slowest, z fastest).
- `modvOriginal`: 1D model table, columns depth (km, positive down), Vp, Vs (km/s).
- `PorS`: 2 = P, 3 = S (uses column `PorS + 1` of the table: Vp or Vs).
- `origin`: grid origin `[lat, lon, altitude]` (degrees, m).

# Returns
- `modv`: n×4 `[x y z v]` (m, km/s).
- `pvel`: velocity as meshgrid array (ny, nx, nz).
- `modvPlot`: n×4 `[lon lat z v]`.
"""
function murat_modv1d(modvXYZ, modvOriginal, PorS, origin)
    modv = hcat(modvXYZ, zeros(size(modvXYZ, 1)))
    # Depth (km, positive down) to altitude (m, positive up), as the grid z
    z1D = -modvOriginal[:, 1] .* 1000
    v1D = modvOriginal[:, PorS+1]
    for k in axes(modv, 1)
        modv[k, 4] = v1D[argmin(abs.(z1D .- modv[k, 3]))]
    end
    xD = sort(unique(modv[:, 1]))
    yD = sort(unique(modv[:, 2]))
    zD = sort(unique(modv[:, 3]), rev = true)
    pvel = fold_field(xD, yD, zD, modv[:, 4])
    modvPlot = plot_coordinates(modv, origin, 0.0)
    return modv, pvel, modvPlot
end

# Lon/lat of grid points in meters from the origin (as MuRAT, with atan(x/y))
"""
    plot_coordinates(modv, origin, az_nan) -> n×4 matrix

Converts the planar coordinates of an n×4 model `[x y z v]` (m from `origin`) to
`[lon lat z v]` (degrees), moving from the origin by distance `sqrt(x²+y²)` and
azimuth `atan(x/y)` on the WGS84 ellipsoid.
MATLAB: the "Change to lat/lon for plotting" blocks of `Murat_modv1D.m` and
`Murat_modv3D.m`, reproduced exactly: `atan(x/y)` (not `atan(x, y)`) is used, so
azimuths are only in (-90°, 90°).

# Arguments
- `modv`: n×4 `[x y z v]`, x east and y north of the origin (m).
- `origin`: `[lat, lon, ...]` (degrees).
- `az_nan`: azimuth used where `x/y` is NaN (node at the origin): 0.0 as in
  `Murat_modv1D`, 45.0 as in `Murat_modv3D` (which also sets the distance to 1 m).
"""
function plot_coordinates(modv, origin, az_nan)
    out = similar(modv)
    for k in axes(modv, 1)
        d = sqrt(modv[k, 1]^2 + modv[k, 2]^2)
        az = rad2deg(atan(modv[k, 1] / modv[k, 2]))
        if isnan(az)
            az = az_nan
            az_nan == 45.0 && (d = 1.0)
        end
        lat2, lon2 = georeckon(origin[1], origin[2], d, az)
        out[k, :] .= (lon2, lat2, modv[k, 3], modv[k, 4])
    end
    return out
end

"""
    murat_modv3d(modvOriginal, origin, ending, nLat, nLon, nzc)

`Murat_modv3D`: 3D model (columns: lat, lon, altitude in m, velocity)
resampled in meters on the inversion grid and on the propagation grid (half
grid step). Returns `(modvP, pvel, modvEqS, modvPlot)`: propagation model as
list and meshgrid, inversion model as list, inversion model in lon/lat.
MATLAB: `Murat_modv3D.m`. Differences: the unused `FPath`, `FLabel`, `mLat`, `mLon`
inputs are dropped and the velocity test figure (`Murat_test_vel_models`) is not made.

# Arguments
- `modvOriginal`: n×4 table `[lat lon altitude v]` (degrees, m, km/s).
- `origin`, `ending`: grid corners `[lat, lon, altitude]` (degrees, m).
- `nLat`, `nLon`, `nzc`: number of inversion nodes in latitude, longitude, vertical.

# Returns
- `modvP`: propagation model, list `[x y z v]` (m, km/s).
- `pvel`: propagation model as meshgrid array (ny, nx, nz).
- `modvEqS`: inversion model, list `[x y z v]` (nLon·nLat·nzc rows).
- `modvPlot`: `modvEqS` with `[lon lat z v]`.
"""
function murat_modv3d(modvOriginal, origin, ending, nLat, nLon, nzc)
    n = size(modvOriginal, 1)
    xyz = zeros(n, 4)
    for k in 1:n
        x, y = local_xy(origin[1], origin[2], modvOriginal[k, 1], modvOriginal[k, 2])
        xyz[k, :] .= (x, y, modvOriginal[k, 3], modvOriginal[k, 4])
    end
    # Original nodes in meters, sorted by x, y and decreasing z (MATLAB sortrows [1,2,-3])
    modvOrC, _ = sortrows_matlab(xyz, [1, 2, -3])
    xOrC, yOrC, zOrC, vOrC = (modvOrC[:, c] for c in 1:4)

    # Inversion grid in meters (z descending), spanning origin to ending
    endX, endY = local_xy(origin[1], origin[2], ending[1], ending[2])
    xM1 = collect(range(min(0, endX), max(0, endX), length = nLon))
    yM1 = collect(range(min(0, endY), max(0, endY), length = nLat))
    zM1 = sort(collect(range(origin[3], ending[3], length = nzc)), rev = true)

    _, _, _, modvEqSpace, _ = rescale_field(xOrC, yOrC, zOrC, vOrC, xM1, yM1, zM1)
    modvPlot = plot_coordinates(modvEqSpace, origin, 45.0)

    # Propagation grid: half the inversion step, built with MATLAB colon semantics
    # so that the number of nodes matches MATLAB exactly
    r2x = abs(xM1[2] - xM1[1]) / 2
    r2y = abs(yM1[2] - yM1[1]) / 2
    r2z = abs(zM1[2] - zM1[1]) / 2
    xp = matlab_colon(xM1[1], r2x, xM1[end])
    yp = matlab_colon(yM1[1], r2y, yM1[end])
    zp = matlab_colon(zM1[1], -r2z, zM1[end])
    _, _, _, modvProp, modvPmesh = rescale_field(xOrC, yOrC, zOrC, vOrC, xp, yp, zp)
    return modvProp, modvPmesh, modvEqSpace, modvPlot
end

# Numeric table separated by spaces, tabs or commas (like MATLAB `load`).
"""
    read_model_table(path) -> Matrix{Float64}

Reads a velocity model file as a numeric matrix (one row per line, with Unix,
Windows or old Mac line endings; separators: spaces, tabs, commas or semicolons;
blank lines skipped).
MATLAB: `load(velocityModel)` in `Murat_checks.m`. Header or comment lines are not
supported (they make `parse` fail); rows of different length raise an error.
"""
function read_model_table(path::AbstractString)
    rows = Vector{Float64}[]
    # Split on \r\n, \r or \n: e.g. iasp91.txt uses carriage returns only
    for line in split(read(path, String), r"\r\n|\r|\n")
        f = split(strip(line), r"[\s,;]+"; keepempty = false)
        isempty(f) && continue
        push!(rows, parse.(Float64, f))
    end
    allequal(length.(rows)) || error("$path: rows have different numbers of columns")
    return permutedims(reduce(hcat, rows))
end
