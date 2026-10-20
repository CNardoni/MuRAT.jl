# Summary values, ParaView files and saving of runs.

"""
    murat_summary(m::Murat)

Values computed by `Murat_plot`: weighted average inverse Qc and its
standard deviation per frequency (`m.Qc.averageQcFrequency`) and the
power-law fit Qc⁻¹ = Q0 f^n (`m.Qc.frequencyFit`).
MATLAB: the `avQcFreq` loop of `Murat_plot.m` and the fit in
`Murat_imageQcFrequency.m` (only displayed in the figure title there), without figures.

# Arguments
- `m`: run after `murat_inversion` (uses `m.Qc.retainQc`, `inverseQc`, `uncertaintyQc`).

# Returns
- `m`, with `m.Qc.averageQcFrequency` (2 × nFreq: weighted mean of retained Qc⁻¹,
  standard deviation) and, for two or more frequencies, `m.Qc.frequencyFit`
  (NamedTuple `Q0`, `Q0_uncertainty`, `exponent`, `exponent_uncertainty`).
"""
function murat_summary(m::Murat)
    cf = m.input.centralFrequency
    nF = length(cf)
    av = zeros(2, nF)
    for k in 1:nF
        # Mean of retained Qc⁻¹ weighted by uncertaintyQc (as MATLAB RZZ weights),
        # plus the unweighted standard deviation
        r = m.Qc.retainQc[:, k]
        Qm = m.Qc.inverseQc[r, k]
        RZZ = m.Qc.uncertaintyQc[r, k]
        av[1, k] = sum(RZZ .* Qm) / sum(RZZ)
        av[2, k] = std(Qm)
    end
    m.Qc.averageQcFrequency = av
    # Linear fit log(Qc⁻¹) = log(Q0) + n log(f) (unweighted, as MATLAB lsqlin(G,d));
    # 95% intervals (1.96σ) from the covariance weighted by 1/std². The Q0
    # uncertainty exp(dm1)*dm1 reproduces the MATLAB figure-title formula.
    if nF >= 2
        G = hcat(ones(nF), log.(cf))
        p = G \ log.(av[1, :])
        covG = inv(G' * Diagonal(av[2, :] .^ -2) * G)
        dm = 1.96 .* sqrt.(diag(covG))
        m.Qc.frequencyFit = (Q0 = exp(p[1]), Q0_uncertainty = exp(dm[1]) * dm[1],
                             exponent = p[2], exponent_uncertainty = dm[2])
    end
    return m
end

# List (MuRAT order, x slowest, z fastest) to an (nx, ny, nz) array
"""
    to_xyz_array(v, nx, ny, nz) -> Array (nx, ny, nz)

Reshapes a field in MuRAT list order (x slowest, z fastest) into an `[i, j, k]`
array indexed by x, y, z, as needed by WriteVTK. No MATLAB counterpart (Julia helper).
"""
to_xyz_array(v, nx, ny, nz) = permutedims(reshape(v, nz, ny, nx), (3, 2, 1))

"""
    murat_paraview(m::Murat)

Writes one ParaView structured grid (`.vts`) per frequency in
`<label>/Results`, with all models and synthetic tests (coordinates: lon,
lat, altitude in km).
No MATLAB counterpart (Julia addition; MATLAB saves `.mat` files of the maps in
`Murat_plot.m`).

# Arguments
- `m`: run after `murat_inversion` (uses `m.PD.modvPeakDelay`, `m.Qc.modvQc`,
  `m.Q.modvQ`, all n×ncol×nFreq).

# Returns
- vector of the written file paths (`<label>/Results/MuRAT_<f>_Hz.vts`, `.` in the
  frequency replaced by `_`).
"""
function murat_paraview(m::Murat)
    g = m.grid
    x, y, z = g.x, g.y, g.z ./ 1000
    nx, ny, nz = length(x), length(y), length(z)
    # Rectilinear node coordinates as 3D arrays (nx, ny, nz): lon, lat, altitude km
    X = [x[i] for i in 1:nx, j in 1:ny, k in 1:nz]
    Y = [y[j] for i in 1:nx, j in 1:ny, k in 1:nz]
    Z = [z[k] for i in 1:nx, j in 1:ny, k in 1:nz]
    files = String[]
    for (k, f) in enumerate(m.input.centralFrequency)
        name = outdir(m, "Results", "MuRAT_$(replace(string(fmt_g(f)), "." => "_"))_Hz")
        vtk_grid(name, X, Y, Z) do vtk
            pd, qc, q = m.PD.modvPeakDelay[:, :, k], m.Qc.modvQc[:, :, k], m.Q.modvQ[:, :, k]
            # Closure: list column -> (nx, ny, nz) array
            A(v) = to_xyz_array(v, nx, ny, nz)
            # Columns of the modv tables: 4 = model value, 5 = hit count,
            # 6/7 = checkerboard input/output, 8/9 = spike input/output
            vtk["PeakDelay"] = A(pd[:, 4])
            vtk["PeakDelayHits"] = A(pd[:, 5])
            vtk["Qc"] = A(qc[:, 4])
            vtk["Q"] = A(q[:, 4])
            for (c, n) in ((6, "CheckIn"), (7, "CheckOut"), (8, "SpikeIn"), (9, "SpikeOut"))
                vtk["Qc_$n"] = A(qc[:, c])
                vtk["Q_$n"] = A(q[:, c])
            end
        end
        push!(files, name * ".vts")
    end
    return files
end

"""
    murat_save(path, m::Murat) -> path

Saves a run with JLD2 (MATLAB: `temp/*.mat`).
The whole `Murat` container is stored under the key `"murat"`, compressed with
Zstandard (lossless; about 6 times smaller than uncompressed for large data sets,
where the inversion matrices take most of the space).
MATLAB: the `save(... 'Murat')` calls of `MuRAT.m`.
"""
murat_save(path::AbstractString, m::Murat) =
    (jldsave(path, JLD2.ZstdFilter(); murat = m); path)

# Fields added to the run by `murat_inversion` and `murat_summary`
const INVERSION_FIELDS = ((:grid, :DiffConstant), (:PD, :modvPeakDelay), (:Qc, :modvQc),
                          (:Qc, :residualRedQc), (:Qc, :dampingValueQc), (:Qc, :lcurves),
                          (:Q, :residualRedQ), (:Q, :dampingValueQ), (:Q, :modvQ),
                          (:Q, :codaConstants), (:Q, :averageQ), (:Q, :estimatedDiffusion),
                          (:Qc, :averageQcFrequency), (:Qc, :frequencyFit))

"""
    murat_save_inversion(path, m::Murat) -> path

Saves only what the inversion stage adds to a run (`INVERSION_FIELDS`: models,
misfit reductions, damping values, L-curves, coda constants, Qc-frequency fit)
and the input, compressed with Zstandard. The data of the forward stage are not
repeated: [`murat_load`](@ref) takes them from `Murat_forward.jld2` in the same
folder. This is how `run_murat` writes `Murat_inverse.jld2` (a few MB).
MATLAB: `save('temp/Murat_inverse.mat','Murat')`, which repeats the whole structure.
"""
function murat_save_inversion(path::AbstractString, m::Murat)
    results = Dict{String,Any}()
    for (sec, f) in INVERSION_FIELDS
        part = getfield(m, sec)
        haskey(part, f) && (results["$sec.$f"] = getproperty(part, f))
    end
    jldsave(path, JLD2.ZstdFilter(); inversion = results, input = m.input)
    return path
end

"""
    murat_load(path) -> Murat

Loads a run saved by `run_murat`. For `Murat_inverse.jld2` (results of the
inversion only) the forward run is read from `Murat_forward.jld2` in the same
folder and the input and inversion results are added to it; other files (and
complete runs saved by older versions or with `murat_save`) are read directly.
Runs saved by an older version of the package also load when `MuratInput` has
gained new fields since: the new fields take their default values (fields that
no longer exist are dropped).
MATLAB: `load('temp/Murat_*.mat')`. The MuRAT package must be loaded (`using MuRAT`)
so that JLD2 can rebuild the `Murat` type.
"""
function murat_load(path::AbstractString)
    upgrade = Dict("MuRAT.MuratInput" => JLD2.Upgrade(MuratInput))
    return jldopen(path, "r"; typemap = upgrade) do f
        haskey(f, "murat") && return f["murat"]
        # Inversion results only: start from the forward run next to this file
        fwd = joinpath(dirname(path), "Murat_forward.jld2")
        isfile(fwd) || error("$path holds only inversion results; $fwd is needed to load the run")
        m = murat_load(fwd)
        inp = f["input"]
        for k in fieldnames(MuratInput)
            setfield!(m.input, k, getfield(inp, k))
        end
        for (key, v) in f["inversion"]
            sec, fld = Symbol.(split(key, "."))
            setproperty!(getfield(m, sec), fld, v)
        end
        return m
    end
end

# Rebuild a stored MuratInput from its fields with the keyword constructor, so that
# fields added after the run was saved get their defaults
function JLD2.rconvert(::Type{MuratInput}, nt::NamedTuple)
    names = fieldnames(MuratInput)
    return MuratInput(; (k => v for (k, v) in pairs(nt) if k in names)...)
end
