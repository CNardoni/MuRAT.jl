# Data selection (Murat_selection, Murat_dataWarning, Murat_components,
# Murat_retainPeakDelay, Murat_retainQc, Murat_declustering).

"""
    murat_data_warning(files, tresholdnoise, maPD, miPD, fT, peakd, Qm, RZZ,
                       rapspcn, rapsp, comp, flag, QcM)

`Murat_dataWarning`: flags data outside the limits and prints the
percentages. Returns `(problemPD, problemQc, problemRZZ, problemQ, yesPD,
compMissing, flag)`.
Called twice by `murat_selection`: on all traces, then on the recordings after
component averaging (the `flag` controls which messages are printed).
MATLAB: `Murat_dataWarning.m`.

# Arguments
- `files`: SAC file names (one per trace).
- `tresholdnoise`: minimum energy ratio (direct/coda and coda/noise).
- `maPD`, `miPD`: maximum and minimum peak delay (s).
- `fT`: threshold on the Qc uncertainty (input `fitTresholdLinear`).
- `peakd`, `Qm`, `RZZ`, `rapspcn`, `rapsp`: peak delay (s), inverse Qc, Qc
  uncertainty, coda-to-noise and direct-to-coda ratios (nData × nFreq each).
- `comp`: number of components; `flag`: 0 on the first call, then the value
  returned by the previous call; `QcM`: `"Linearized"` or `"NonLinear"`.

# Returns
- `problemPD`, `problemQc`, `problemRZZ`, `problemQ`: per frequency, the file
  names failing each test (vector of nFreq vectors).
- `yesPD`: peak delays inside (miPD, maPD) (nData × nFreq Bool).
- `compMissing`: nData × nFreq × 3 Bool, failing peak delay / Qc / energy ratios
  (used by `murat_components` to skip those components).
- `flag`: updated message flag (2 after a one-component call, 1 after the first
  multi-component call).
"""
function murat_data_warning(files, tresholdnoise, maPD, miPD, fT, peakd, Qm, RZZ,
                            rapspcn, rapsp, comp, flag, QcM)
    nRows, nCols = size(Qm)
    # Bool masks (nData × nFreq) of the data failing each test; the Qc uncertainty
    # test has opposite directions for the two Qc methods, as in MATLAB
    yes_pd = (peakd .> miPD) .& (peakd .< maPD)
    no_pd = .!yes_pd
    no_Qc = Qm .== 0
    no_RZZ = QcM == "Linearized" ? RZZ .<= fT : RZZ .< fT
    no_Q = (rapspcn .< tresholdnoise) .| (rapsp .< tresholdnoise)
    compMissing = cat(no_pd, no_Qc, no_Q, dims = 3)
    # findall: like MATLAB files(mask), also when mask is shorter than files
    # (second call, after averaging the components)
    pick(mask) = [files[findall(mask[:, c])] for c in 1:nCols]
    # Percentage of failing data per frequency
    pct(mask) = vec(sum(mask, dims = 1)) ./ nRows .* 100
    pQc, pRZZ, pPD, pQ = pct(no_Qc), pct(no_RZZ), pct(no_pd), pct(no_Q)
    # One component: print per-frequency percentages once (flag -> 2). More
    # components: print the means over frequencies on both calls
    if comp == 1 && flag != 2
        println("$pQc % of Qc are below threshold")
        println("$pRZZ % of Qc uncertainties are below threshold")
        println("$pPD % of peak delays are outside limits")
        println("$pQ % of ratios are below threshold")
        flag = 2
    elseif flag in (0, 1)
        println("$(mean(pQc)) % of Qc are below threshold")
        println("$(mean(pRZZ)) % of Qc uncertainties are below threshold")
        println("$(mean(pPD)) % of peak delays are outside limits")
        println("$(mean(pQ)) % of ratios are below threshold")
        if flag == 0
            println("Processing to see how many data you have when considering $comp components")
            flag = 1
        end
    end
    return pick(no_pd), pick(no_Qc), pick(no_RZZ), pick(no_Q), yes_pd, compMissing, flag
end

"""
    murat_components(comp, peakd, Qm, RZZ, rapsp, rapspcn, compMissing)

`Murat_components`: averages the measurements of the components of each
recording, ignoring the missing ones.
The components of one recording are consecutive traces (`comp` rows per
group); a value is NaN when no component of the group is valid.
MATLAB: `Murat_components.m`.

# Arguments
- `comp`: number of components (2 or 3).
- `peakd`, `Qm`, `RZZ`, `rapsp`, `rapspcn`: per-trace values (nData × nFreq).
- `compMissing`: nData × nFreq × 3 Bool from `murat_data_warning` (peak delay,
  Qc, energy ratios); Qc and its uncertainty share the Qc mask.

# Returns
- `(peakd, Qm, RZZ, rapsp, rapspcn)` averaged per recording ((nData ÷ comp) × nFreq).
"""
function murat_components(comp, peakd1, Qm1, RZZ1, rapsp1, rapspcn1, compMissing)
    dataL, lcf = size(peakd1)
    nG = dataL ÷ comp
    peakd, Qm, RZZ, rapsp, rapspcn = (fill(NaN, nG, lcf) for _ in 1:5)
    # Mean of the valid entries of each column of a comp × nFreq block, written to
    # row g of `out` (NaN where no entry is valid)
    function avg!(out, blk, valid, g)
        n = vec(sum(valid, dims = 1))
        s = vec(sum(blk .* valid, dims = 1))
        out[g, :] .= ifelse.(n .== 0, NaN, s ./ n)
    end
    for (g, i) in enumerate(1:comp:dataL)
        rows = i:i+comp-1
        vPD = .!compMissing[rows, :, 1]
        vQc = .!compMissing[rows, :, 2]
        vR = .!compMissing[rows, :, 3]
        avg!(peakd, peakd1[rows, :], vPD, g)
        avg!(Qm, Qm1[rows, :], vQc, g)
        avg!(RZZ, RZZ1[rows, :], vQc, g)
        avg!(rapsp, rapsp1[rows, :], vR, g)
        avg!(rapspcn, rapspcn1[rows, :], vR, g)
    end
    return peakd, Qm, RZZ, rapsp, rapspcn
end

"""
    murat_retain_peak_delay(t_phase, l10pd, yesPD, Apd, stDevPD)

`Murat_retainPeakDelay`: robust line fit of log10(peak delay) against
log10(travel time) (all recordings with a valid peak delay), residuals and outlier
removal at `stDevPD` standard deviations.
The residuals (peak-delay variations) are the data of the peak-delay
(scattering) inversion. Called once per frequency.
MATLAB: `Murat_retainPeakDelay.m` (`fit(...,'poly1','Robust','on')` is
reproduced by `robust_line_fit`). MuRAT 4.0 fits all recordings with a valid peak
delay; MuRAT 3 also excluded travel times below 2.6 s.

# Arguments
- `t_phase`: log10 of the travel time (s), one per recording.
- `l10pd`: log10 of the peak delay (s), one per recording.
- `yesPD`: Bool, peak delay inside the allowed limits.
- `Apd`: peak-delay inversion matrix (nRecordings × nModel).
- `stDevPD`: number of standard deviations for the outlier bounds.

# Returns
- `pab`: `[slope, intercept]` of the robust line.
- `lpdelta`: residuals log10(pd) - (slope t_phase + intercept), one per recording.
- `retain`: Bool, recordings used in the inversion.
- `rayCrosses`: Bool per cell, crossed by at least one retained ray.
"""
function murat_retain_peak_delay(t_phase, l10pd, yesPD, Apd, stDevPD)
    yes = yesPD
    pab = robust_line_fit(t_phase[yes], l10pd[yes])
    lpdelta = l10pd .- (pab[1] .* t_phase .+ pab[2])
    # Mean and standard deviation of the residuals ignoring NaN (MATLAB 'omitnan')
    vals = filter(!isnan, lpdelta[yes])
    mu, sigma = mean(vals), std(vals)
    lB, uB = mu - stDevPD * sigma, mu + stDevPD * sigma
    # Keep the selected recordings whose residual lies inside [lB, uB]
    retain = falses(length(yes))
    sel = findall(yes)
    retain[sel] .= .!((lpdelta[sel] .< lB) .| (lpdelta[sel] .> uB))
    rayCrosses = vec(sum(Apd[retain, :], dims = 1)) .!= 0
    return pab, lpdelta, retain, rayCrosses
end

"""
    murat_retain_qc(fT, Qm, RZZ, QcM) -> retain

`Murat_retainQc`: Qc measurements kept for the inversion.
`"Linearized"`: positive Qc, `RZZ < 1/fT` (i.e. the correlation coefficient |R| of
the fit above `fT`), and Qc below mean + 2 standard deviations of these.
`"NonLinear"`: positive Qc and uncertainty below `fT` (below its maximum when
`fT = 0`). Called once per frequency.
MATLAB: `Murat_retainQc.m` of MuRAT 4.0 (MuRAT 3 also required rays shorter than
100 km and used mean + 1 standard deviation).

# Arguments
- `fT`: uncertainty threshold (input `fitTresholdLinear`).
- `Qm`, `RZZ`: inverse Qc and its uncertainty, one per recording.
- `QcM`: `"Linearized"` or `"NonLinear"`.

# Returns
- `retain`: Bool, one per recording.
"""
function murat_retain_qc(fT, Qm, RZZ, QcM)
    if QcM == "Linearized"
        tmp = (Qm .> 0) .& (RZZ .< 1 / fT)
        lim = mean(Qm[tmp]) + 2 * std(Qm[tmp])
        return tmp .& (Qm .< lim)
    elseif QcM == "NonLinear"
        # fT = 0: below the largest uncertainty (MATLAB max ignores NaN)
        return (Qm .> 0) .& (RZZ .< (fT == 0 ? maximum(filter(!isnan, RZZ)) : fT))
    end
    error("Unknown Qc measurement $QcM")
end

"""
    murat_declustering(origin, ending, x, y, z, QcM, uncertainty, locationsDeg, factor)

`Murat_declustering`: divides the grid in cells `factor` times smaller and,
in each cell, keeps for each station the event-station pair with the
smallest Qc uncertainty.
Used to reduce the weight of event clusters, for Qc (with the Qc uncertainty)
and for Q (with 1/coda-to-noise ratio as uncertainty). Called once per
frequency when the input `declustering` is set.
MATLAB: `Murat_declustering.m`. The MATLAB quirks are kept: the cell is located
from the event position, stations are identified by their latitude only, the
`"Linear"` test never matches `"Linearized"` (so 1/uncertainty is always used),
and the latitude/longitude step counts follow the original code (see README).

# Arguments
- `origin`, `ending`: grid corners `[lat, lon, altitude (m)]`.
- `x`, `y`, `z`: grid node vectors (only their lengths are used).
- `QcM`: Qc method name; `uncertainty`: one value per recording (smaller is better).
- `locationsDeg`: `[evla evlo evdp(m) stla stlo stel]`, one row per recording.
- `factor`: refinement factor of the grid (input `declustering`).

# Returns
- `indKeep`: Bool, one per recording, the pairs kept.
"""
function murat_declustering(origin, ending, x, y, z, QcM, uncertainty, locationsDeg, factor)
    # Refined grid. lon_grid spans origin[1] -> ending[1] and is compared with column 1
    # (evla), lat_grid with column 2 (evlo): names swapped, values consistent, as in MATLAB.
    # z_grid is depth (m, positive down), like evdp in column 3
    lat_step = length(y) * factor
    lon_step = length(x) * factor
    z_step = length(z) * factor
    lat_grid = collect(range(origin[2], ending[2], length = lat_step))
    lon_grid = collect(range(origin[1], ending[1], length = lon_step))
    z_grid = -collect(range(origin[3], ending[3], length = z_step))
    # One row per event-station pair; f below indexes these unique rows
    loc = unique_rows_stable(locationsDeg)
    indKeep = falses(length(uncertainty))
    kept = zeros(0, 8)
    # For each refined cell: events strictly inside it; append 1/uncertainty (column 7)
    # and row index (column 8), sort by station latitude then column 7, and keep the
    # last row per station, i.e. the largest 1/uncertainty
    for i in 1:lon_step-1, j in 1:lat_step-1, k in 1:z_step-1
        f = findall((loc[:, 1] .> lon_grid[i]) .& (loc[:, 1] .< lon_grid[i+1]) .&
                    (loc[:, 2] .> lat_grid[j]) .& (loc[:, 2] .< lat_grid[j+1]) .&
                    (loc[:, 3] .> z_grid[k]) .& (loc[:, 3] .< z_grid[k+1]))
        isempty(f) && continue
        ev = hcat(loc[f, :], QcM == "Linear" ? uncertainty[f] : 1 ./ uncertainty[f], f)
        d, _ = sortrows_matlab(ev, [4, 7])
        ia = unique_last_indices(d[:, 4])
        kept = vcat(kept, d[ia, :])
    end
    isempty(kept) && return indKeep
    kept, _ = sortrows_matlab(kept, [8])
    indKeep[Int.(kept[:, 8])] .= true
    return indKeep
end

"""
    murat_selection(m::Murat)

`Murat_selection`: warnings, component averaging, peak-delay residuals,
retained data and crossed cells for the three inversions, declustering.
Called at the end of `murat_data`; prepares the data and masks used by the
peak-delay, Qc and Q inversions.
MATLAB: `Murat_selection.m` (the declustering figures are not produced here;
the declustering masks are stored in `m.Qc.declustered` and `m.Q.declustered`).

# Arguments
- `m`: `Murat` container filled by `murat_data` (one row per trace).

# Returns
- `m`, updated in place with one row per recording (nRec = nData / components):
  averaged measurements, travel time (s), `m.PD.fitrobust` (2 × nFreq),
  `m.PD.variationPeakDelay` (nRec × nFreq), the retained-data masks
  (`retainPeakDelay`, `retainQc`, `retainQ`, nRec × nFreq), the crossed-cell
  masks (`raysPeakDelay`, `raysQc`, `raysQ`, nModel × nFreq), the reduced
  inversion matrices and the file names of the rejected data.
"""
function murat_selection(m::Murat)
    s = m.input
    g = m.grid
    components = s.components
    QcM = s.QcMeasurement
    cf = s.centralFrequency
    fT = s.fitTresholdLinear
    Declust = s.declustering

    Apd = m.PD.inversionMatrixPeakDelay
    A = m.Q.inversionMatrixQ
    Ac = m.Qc.inversionMatrixQc
    peakd = m.PD.peakDelay
    luntot = m.rays.totalLengthRay
    tPS = m.rays.theoreticalTime
    locationsDeg = m.rays.locationsDeg
    Qm = m.Qc.inverseQc
    RZZ = m.Qc.uncertaintyQc
    rapsp = m.Q.energyRatioBodyCoda
    rapspcn = m.Q.energyRatioCodaNoise
    raysplot = m.rays.raysPlot
    tCoda = m.Qc.tCoda
    dataL, nFreq = size(peakd)
    nModel = size(g.modv, 1)

    Declust === nothing || println("Waveforms used before Qc and Q declustering: $(length(g.listSac))")

    # First warning on all traces: problem lists and missing components
    _, problemQc, problemRZZ, problemQ, _, compMissing, flag =
        murat_data_warning(g.listSac, s.tresholdNoise, s.maximumPeakDelay, s.minimumPeakDelay,
                           fT, peakd, Qm, RZZ, rapspcn, rapsp, components, 0, QcM)

    # Keep one row per recording (first component; rays and kernels were computed there)
    idx = 1:components:dataL
    luntot = luntot[idx]
    time0 = tPS[idx]
    tCoda = tCoda[idx, :]
    locationsDeg = locationsDeg[idx, :]
    raysplot = raysplot[:, :, idx]
    Ac = Ac[idx, :, :]
    Apd = Apd[idx, :]
    A = A[idx, :]
    if components > 1
        peakd, Qm, RZZ, rapsp, rapspcn = murat_components(components, peakd, Qm, RZZ,
                                                          rapsp, rapspcn, compMissing)
    end
    # Second warning on the recordings: valid peak delays and the peak-delay problem list
    problemPD, _, _, _, yesPD, _, _ =
        murat_data_warning(g.listSac, s.tresholdNoise, s.maximumPeakDelay, s.minimumPeakDelay,
                           fT, peakd, Qm, RZZ, rapspcn, rapsp, components, flag, QcM)

    # With P picks the travel time is scaled by Vp/Vs = sqrt(3) to approximate the S
    # travel time (peak delays are measured on the S-wave maximum)
    t_phase = s.POrS == 2 ? log10.(time0 .* sqrt(3)) : log10.(time0)
    s.POrS == 2 && (time0 = time0 .* sqrt(3))
    l10pd = log10.(peakd)

    n = size(peakd, 1)
    fitrobust = zeros(2, nFreq)
    lpdelta = zeros(n, nFreq)
    retain_pd = falses(n, nFreq)
    retain_Qm = falses(n, nFreq)
    retain_Q = falses(n, nFreq)
    ray_pd = falses(nModel, nFreq)
    ray_Qc = falses(nModel, nFreq)
    ray_Q = falses(nModel, nFreq)
    declusteredQc = Vector{Any}(undef, nFreq)
    declusteredQ = Vector{Any}(undef, nFreq)

    # Per frequency: peak-delay residuals and outliers, Qc and Q retained data
    # (optionally declustered) and the cells crossed by the retained data
    for i in 1:nFreq
        pab, lpd, rpd, rcpd = murat_retain_peak_delay(t_phase, l10pd[:, i], yesPD[:, i], Apd, s.stDevPD)
        rQc = murat_retain_qc(fT, Qm[:, i], RZZ[:, i], QcM)
        if Declust !== nothing
            keep = murat_declustering(s.origin, s.ending, g.x, g.y, g.z, QcM, RZZ[:, i],
                                      locationsDeg, Declust)
            rQc = rQc .& keep
            declusteredQc[i] = keep
        end
        # A cell counts for Qc if the summed kernel weight of retained data exceeds 1e-19
        rcQc = vec(sum(Ac[rQc, :, i], dims = 1)) .> 1e-19

        # Q (coda normalization): valid Qc and both energy ratios above the noise threshold
        rQ = .!isnan.(Qm[:, i]) .& (rapspcn[:, i] .>= s.tresholdNoise) .& (rapsp[:, i] .>= s.tresholdNoise)
        if Declust !== nothing
            keep = murat_declustering(s.origin, s.ending, g.x, g.y, g.z, QcM, 1 ./ rapspcn[:, i],
                                      locationsDeg, Declust)
            rQ = rQ .& keep
            declusteredQ[i] = keep
        end
        rcQ = vec(sum(A[rQ, :], dims = 1)) .!= 0

        fitrobust[:, i] .= pab
        lpdelta[:, i] .= lpd
        retain_pd[:, i] .= rpd
        ray_pd[:, i] .= rcpd
        retain_Qm[:, i] .= rQc
        ray_Qc[:, i] .= rcQc
        retain_Q[:, i] .= rQ
        ray_Q[:, i] .= rcQ
    end
    println("Waveforms used after Qc declustering: ", vec(sum(retain_Qm, dims = 1)))
    println("Waveforms used after Q declustering: ", vec(sum(retain_Q, dims = 1)))

    m.rays.locationsDeg = locationsDeg
    m.rays.totalLengthRay = luntot
    m.rays.raysPlot = raysplot
    m.rays.travelTime = time0
    m.PD.peakDelay = peakd
    m.PD.problemPD = problemPD
    m.PD.variationPeakDelay = lpdelta
    m.PD.fitrobust = fitrobust
    m.PD.retainPeakDelay = retain_pd
    m.PD.raysPeakDelay = ray_pd
    m.PD.inversionMatrixPeakDelay = Apd
    m.Qc.tCoda = tCoda
    m.Qc.inverseQc = Qm
    m.Qc.uncertaintyQc = RZZ
    m.Qc.problemQc = problemQc
    m.Qc.problemRZZ = problemRZZ
    m.Qc.raysQc = ray_Qc
    m.Qc.retainQc = retain_Qm
    m.Qc.inversionMatrixQc = Ac
    Declust === nothing || (m.Qc.declustered = declusteredQc; m.Q.declustered = declusteredQ)
    m.Q.problemQ = problemQ
    m.Q.energyRatioBodyCoda = rapsp
    m.Q.energyRatioCodaNoise = rapspcn
    m.Q.retainQ = retain_Q
    m.Q.raysQ = ray_Q
    m.Q.inversionMatrixQ = A
    return m
end
