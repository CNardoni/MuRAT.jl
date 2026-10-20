# MATLAB `.mat` output of a MuRAT run (package extension, loaded with `using MAT`).
# Adds the method `MuRAT.murat_matfiles`, called by `run_murat` when the input has
# `SaveMatFiles = 1`.
module MuratMATExt

using MuRAT
using MuRAT: Murat, MuratInput, MuratSection, outdir, fmt_g, fold_field, meshgrid3,
             interp3_linear
using MAT

# Arrays with more elements than this are left out of Murat_summary.mat
# (inversion matrices, ray paths, ray crossings: GB for large data sets)
const MAX_ELEMENTS = 2_000_000

"""
    fcname(f) -> String

Frequency label used in the file names, as MATLAB `strrep(num2str(f),'.','_')`.
"""
fcname(f) = replace(string(fmt_g(f)), "." => "_")

"""
    refined_grid(x, y, z, divi) -> (xp, yp, zp_km)

Query vectors of `Murat_image3D.m`: grid steps divided by `divi`, one extra step
beyond each end of the grid, same number of points as MATLAB `linspace`.

# Arguments
- `x`, `y`: grid longitudes and latitudes (°); `z`: grid altitudes (m).
- `divi`: refinement factor (5 in `Murat_plot.m`).

# Returns
- `xp`, `yp` (°) and `zp_km` (km) refined vectors.
"""
function refined_grid(x, y, z, divi)
    dx, dy, dz = (x[2] - x[1]) / divi, (y[2] - y[1]) / divi, (z[2] - z[1]) / divi
    n(a, d) = round(Int, (a[end] - a[1]) / d, RoundNearestTiesAway) + 1
    xp = collect(range(x[1] - dx, x[end] + dx, length = n(x, dx)))
    yp = collect(range(y[1] - dy, y[end] + dy, length = n(y, dy)))
    zp = collect(range(z[1] - dz, z[end] + dz, length = n(z, dz)))
    return xp, yp, zp ./ 1000
end

"""
    interpolated(x, y, z, V, divi) -> (Xi, Yi, Zi, Vi)

Field `V` (meshgrid order, ny × nx × nz) linearly interpolated on the refined grid
of `Murat_image3D.m` (NaN outside the grid), before its Gaussian smoothing: the
`*_inter_*.mat` files of `Murat_plot.m`.

# Arguments
- `x`, `y`: grid longitudes and latitudes (°); `z`: grid altitudes (m).
- `V`: field on the inversion grid; `divi`: refinement factor.

# Returns
- `Xi`, `Yi`, `Zi` (°, °, km) meshgrid arrays and the interpolated field `Vi`.
"""
function interpolated(x, y, z, V, divi)
    xp, yp, zp = refined_grid(x, y, z, divi)
    Xi, Yi, Zi = meshgrid3(xp, yp, zp)
    Vi = interp3_linear(x, y, z ./ 1000, V, Xi, Yi, Zi)
    return Xi, Yi, Zi, Vi
end

"""
    tomat(v) -> value MAT.jl can write, or `nothing` to leave it out

Converts a Julia value of the `Murat` container to a MATLAB-compatible one:
numbers, strings and numeric arrays as they are, `nothing` as `[]`, symbols as
strings, long string lists as char matrices, structs, `Dict`s and `NamedTuple`s as
MATLAB structures. Arrays larger than `MAX_ELEMENTS` and unsupported values return
`nothing`; their names are collected in `skipped`.
"""
function tomat(v, name, skipped)
    v === nothing && return zeros(0, 0)
    v isa Union{Number,AbstractString} && return v
    v isa Symbol && return String(v)
    if v isa AbstractArray
        if length(v) > MAX_ELEMENTS
            push!(skipped, "$name ($(join(size(v), "×")))")
            return nothing
        end
        eltype(v) <: Bool && return Array{Bool}(v)
        eltype(v) <: Number && return Array(v)
        if all(e -> e isa AbstractString, v)
            # Long lists (file names) as a char matrix: one name per row
            length(v) <= 5000 && return String.(collect(v))
            L = maximum(length, v; init = 0)
            return permutedims(reduce(hcat, [collect(rpad(s, L)) for s in v]))
        end
        # Vector of numeric tuples (e.g. groups of SAC header values): one row each
        if v isa AbstractVector && !isempty(v) && all(e -> e isa Tuple && all(x -> x isa Number, e), v)
            return permutedims(reduce(hcat, [collect(e) for e in v]))
        end
        # Vector of structs (e.g. SAC headers): one array per field
        if v isa AbstractVector && !isempty(v) && isstructtype(eltype(v)) &&
           !(eltype(v) <: Union{Tuple,AbstractArray,AbstractString}) && fieldcount(eltype(v)) > 0
            return tomat(Dict(String(f) => [getfield(e, f) for e in v]
                              for f in fieldnames(eltype(v))), name, skipped)
        end
        # Other arrays: MATLAB cell array of converted elements
        c = Any[tomat(e, "$name{$i}", skipped) for (i, e) in enumerate(v)]
        any(isnothing, c) && return nothing
        return reshape(c, size(v))
    end
    if v isa AbstractDict || v isa NamedTuple || v isa MuratSection
        keyz = v isa MuratSection ? propertynames(v) : keys(v)
        d = Dict{String,Any}()
        for k in keyz
            # A field that cannot be converted is left out, not the whole structure
            e = try
                tomat(v isa MuratSection ? getproperty(v, k) : v[k], "$name.$k", skipped)
            catch err
                push!(skipped, "$name.$k ($(typeof(err)))")
                nothing
            end
            e === nothing || (d[string(k)] = e)
        end
        return d
    end
    v isa Tuple && return tomat(collect(Any, v), name, skipped)
    if isstructtype(typeof(v)) && fieldcount(typeof(v)) > 0
        return tomat(Dict(String(f) => getfield(v, f) for f in fieldnames(typeof(v))),
                     name, skipped)
    end
    push!(skipped, "$name ($(typeof(v)))")
    return nothing
end

"""
    murat_summary_struct(m::Murat) -> (Dict, skipped)

The `Murat` container as a MATLAB structure (input, grid, rays, PD, Qc, Q),
without the arrays larger than `MAX_ELEMENTS`.

# Returns
- the structure as nested `Dict`s, and the list of fields left out.
"""
function murat_summary_struct(m::Murat)
    skipped = String[]
    d = Dict{String,Any}()
    for sec in (:input, :grid, :rays, :PD, :Qc, :Q)
        e = tomat(getfield(m, sec), "Murat.$sec", skipped)
        e === nothing || (d[String(sec)] = e)
    end
    # Record what was left out, so MATLAB users know it is not an error
    d["skippedFields"] = isempty(skipped) ? String[] : skipped
    return d, skipped
end

"""
    MuRAT.murat_matfiles(m::Murat)

Writes `<label>/matfiles/*.mat` as `Murat_plot.m` (per frequency: `PD_inter`,
`PD_not_inter`, `Qc_inter`, `Qc_not_inter`, `Q`, `hitmap`, `Qc_check_in`,
`Qc_check_out`, with the MATLAB variable names) and `<label>/Murat_summary.mat`
(variable `Murat`: the run without its largest arrays, listed in
`Murat.skippedFields`).
MATLAB: the `save` calls of `Murat_plot.m`. Differences: the files are written
whatever the `Plot*` flags; `Murat_summary.mat` replaces the full `Murat.mat`.

# Arguments
- `m`: a run after `murat_inversion` and `murat_summary`.
"""
function MuRAT.murat_matfiles(m::Murat)
    s, g = m.input, m.grid
    x, y, z = g.x, g.y, g.z
    folder = outdir(m, "matfiles")
    mkpath(folder)
    save(name, vars...) = matwrite(joinpath(folder, name), Dict{String,Any}(vars...))

    # Meshgrid coordinates of Murat_fold (Z in km, as Murat_plot.m)
    X, Y, Z1 = meshgrid3(x, y, z)
    Z = Z1 ./ 1000
    fold(v) = fold_field(x, y, z, v)
    # Events and stations: lat, lon, depth (km) as in Murat_plot.m
    ev = m.rays.locationsDeg
    evst = hcat(ev[:, 1:2], -ev[:, 3] ./ 1000, ev[:, 4:5], ev[:, 6] ./ 1000)
    divi = 5

    for (k, f) in enumerate(s.centralFrequency)
        fc = fcname(f)
        mPD = fold(m.PD.modvPeakDelay[:, 4, k])
        mQc = fold(m.Qc.modvQc[:, 4, k])
        mQ = fold(m.Q.modvQ[:, 4, k])

        # Peak delays reduced to the cells above 5% of the extreme values (MuRAT 4.0; 3% in MuRAT 3)
        pos, neg = filter(>=(0), mPD), filter(<(0), mPD)
        thrP = isempty(pos) ? 0.0 : maximum(pos) / 100 * 5
        thrN = isempty(neg) ? 0.0 : minimum(neg) / 100 * 5
        mPDRed = mPD .* ((mPD .> thrP) .| (mPD .< thrN))

        Xi, Yi, Zi, pd_inter = interpolated(x, y, z, mPDRed, divi)
        save("PD_inter_$(fc)_Hz.mat", "Xi" => Xi, "Yi" => Yi, "Zi" => Zi, "pd_inter" => pd_inter)
        save("PD_not_inter_$(fc)_Hz.mat", "X" => X, "Y" => Y, "Z" => Z, "mPDRed" => mPDRed)

        Xi, Yi, Zi, qc_inter = interpolated(x, y, z, mQc, divi)
        save("Qc_inter_$(fc)_Hz.mat", "Xi" => Xi, "Yi" => Yi, "Zi" => Zi, "qc_inter" => qc_inter)
        save("Qc_not_inter_$(fc)_Hz.mat", "X" => X, "Y" => Y, "Z" => Z, "mQc" => mQc)

        save("Q_$(fc)_Hz.mat", "X" => X, "Y" => Y, "Z" => Z, "mQ" => mQ)

        # Hit map: grid nodes, peak-delay ray count, retained events/stations
        save("hitmap_$(fc)_Hz.mat", "CoordGrid" => g.DDcoordinates,
             "rPD" => m.PD.modvPeakDelay[:, 5, k],
             "evst_pd" => evst[m.PD.retainPeakDelay[:, k], :])

        # Qc checkerboard input (column 6) and output (column 7)
        save("Qc_check_in$(fc)_Hz.mat", "X" => X, "Y" => Y, "Z" => Z,
             "check_inputQc" => fold(m.Qc.modvQc[:, 6, k]))
        save("Qc_check_out$(fc)_Hz.mat", "X" => X, "Y" => Y, "Z" => Z,
             "check_outputQc" => fold(m.Qc.modvQc[:, 7, k]),
             "evst_Qc" => evst[m.Qc.retainQc[:, k], :], "x" => x, "y" => y, "z" => z)
    end

    # Compact Murat structure
    d, skipped = murat_summary_struct(m)
    matwrite(outdir(m, "Murat_summary.mat"), Dict{String,Any}("Murat" => d))
    isempty(skipped) ||
        @info "Murat_summary.mat: large or unsupported fields left out" skipped
    return nothing
end

end # module
