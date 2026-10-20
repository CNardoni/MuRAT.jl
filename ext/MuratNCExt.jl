# LITHO1.0 velocity models (NetCDF), package extension loaded with `using NCDatasets`.
# Adds the method `MuRAT.read_litho1`, called by `murat_checks` when
# `availableVelocity = 0` and `namev = "LITHO1.0.nc"`.
module MuratNCExt

using MuRAT
using MuRAT: pchip
using NCDatasets
using Statistics

"""
    MuRAT.read_litho1(path, lat, lon) -> Matrix{Float64}

LITHO1.0 profile at the node nearest to (`lat`, `lon`), translated from
`Murat_readLithos1.m`: for every layer boundary `<name>_depth` of the file, the
depth and the boundary properties (`vp`, `vs`, ...) at that node are collected;
for each property the (depth, value) pairs are sorted by depth; Vp and Vs are
averaged at equal depths (a zero is replaced by half the next value), a point is
added at 5 km above sea level, and the profile is interpolated every 1 km with
`pchip`.
MATLAB: `Murat_readLithos1.m`. As in MATLAB the profile values are kept in single
precision (the file is single precision) until the final table.

# Returns
- table `[depth (km), 6371 - depth, vp, vs]` (km, km, km/s, km/s).
"""
function MuRAT.read_litho1(path::String, lat::Float64, lon::Float64)
    NCDataset(path) do ds
        lats, lons = Array(ds["latitude"]), Array(ds["longitude"])
        iLat, iLon = argmin(abs.(lats .- lat)), argmin(abs.(lons .- lon))
        names = keys(ds)
        # Value at the node (lon, lat) of a 2D variable; fill values become NaN
        at(name) = (v = ds[name][iLon, iLat]; ismissing(v) ? NaN32 : Float32(v))
        props = Dict{String,Vector{Tuple{Float32,Float32}}}()
        for dname in filter(n -> endswith(n, "_depth"), names)
            base = dname[1:end-length("_depth")]
            depth = at(dname)
            for vname in names
                (startswith(vname, base * "_") && !endswith(vname, "_depth")) || continue
                prop = vname[length(base)+2:end]
                push!(get!(props, prop, Tuple{Float32,Float32}[]), (depth, at(vname)))
            end
        end
        # Per property: drop NaN depths, stable sort by depth (as MATLAB sort)
        clean(p) = sort(filter(t -> !isnan(t[1]), props[p]), by = first, alg = MergeSort)
        vp, vs = clean("vp"), clean("vs")
        tab = hcat(first.(vp), 6371 .- first.(vp), last.(vp), last.(vs))
        # sortrows by depth then vp; mean of vp and vs at equal depths
        tab = tab[sortperm(collect(eachrow(tab[:, [1, 3]])), by = r -> (r[1], r[2])), :]
        depths = unique(tab[:, 1])
        m3 = [mean(tab[tab[:, 1].==dv, 3]) for dv in depths]
        m4 = [mean(tab[tab[:, 1].==dv, 4]) for dv in depths]
        # A zero (not in the last row) becomes half the next value
        for m in (m3, m4)
            idx = findall(i -> m[i] == 0 && i < length(m), eachindex(m))
            m[idx] .= m[idx.+1] ./ 2
        end
        x = vcat(Float32(-5), depths)
        x3, x4 = vcat(m3[1], m3), vcat(m4[1], m4)
        xq = collect(minimum(x):1:maximum(x))
        y3, y4 = pchip(x, x3, xq), pchip(x, x4, xq)
        return Float64.(hcat(xq, 6371 .- xq, y3, y4))
    end
end

end # module
