# Murat_data / Murat_dataParallelized: measurements and forward modelling
# for every trace.

# Traces that carry the ray-dependent computations (first of each component group)
"""
    recognize_components(i, components) -> Bool

`true` if trace `i` must carry the ray-dependent computations (ray tracing,
peak-delay/Q matrix rows, Qc kernel). With 2 or 3 components the files of one
recording are consecutive, so only the first of each group (i = 1, 1 + comp,
...) is processed; with 1 component every trace is; any other value disables them.
MATLAB: local function `recognizeComponents` in `Murat_dataParallelized.m`.
"""
recognize_components(i, components) =
    components == 1 ? true : (components in (2, 3) ? mod(i, components) == 1 : false)

"""
    murat_data(m::Murat) -> Murat

Measures peak delay, Qc and energy ratios on each trace, traces the rays and
computes the kernels (`Murat_dataParallelized`). Uses all Julia threads
(start Julia with `-t N`). Then runs the data selection (`murat_selection`).
This is the data step of the workflow, after `murat_checks` (grid, velocity
model, SAC list and headers) and before the inversions.
MATLAB: `Murat_dataParallelized.m` (and the serial `Murat_data.m`). Unlike
MATLAB, no kernel figures are saved here; the kernel of the first trace is
stored in `m.Qc.firstKernel` instead.

The Qc kernels follow `kernelType`: `"differential"` (MuRAT 4.0, default) uses the
kernels at the start and at the end of the coda window (`murat_coda_matrix_differential`),
`"middle"` (MuRAT 3) one kernel at the middle of the window (`murat_coda_matrix`).

Traces that fail the coda-length or travel-time checks keep NaN measurements
(and zero Qc matrix rows); their count is printed as a percentage at the end.

# Arguments
- `m`: `Murat` container after `murat_checks` (uses `m.input` and `m.grid`).

# Returns
- `m`, filled in place. Per trace (nData = number of SAC files, nFreq = number
  of central frequencies, nModel = grid cells): `m.rays` (locations in degrees
  and m, nData × 6; travel time pick - origin, s; ray length, m;
  ray points for plotting, 100 × 5 × nData; hit count per cell), `m.PD`
  (peak delay, s, nData × nFreq; matrix nData × nModel), `m.Q` (energy
  ratios, nData × nFreq; matrix nData × nModel), `m.Qc` (inverse Qc and
  uncertainty, nData × nFreq; lapse time of coda start, s; matrix
  nData × nModel × nFreq), `m.grid.modvQc`. `murat_selection` then reduces
  these to one row per recording and adds the retained-data masks.
"""
function murat_data(m::Murat)
    s = m.input
    g = m.grid
    files = g.listSac
    nData = length(files)
    cf = s.centralFrequency
    nCF = length(cf)
    modv = g.modv
    nModel = size(modv, 1)
    # Qc kernels are evaluated at the cell centres: grid nodes shifted by half a step
    # (modv columns: x, y, z in m)
    gridStepZ = (modv[2, 3] - modv[1, 3]) / 2
    step = (g.gridStepX / 2, g.gridStepY / 2, gridStepZ)
    modvQc = hcat(modv[:, 1] .+ step[1], modv[:, 2] .+ step[2], modv[:, 3] .+ step[3])

    # Outputs: one row per trace; NaN stays where the trace is rejected
    locationDeg = zeros(nData, 6)
    locationM = zeros(nData, 6)
    theoreticalTime = zeros(nData)
    totalLengthRay = zeros(nData)
    peakDelay = fill(NaN, nData, nCF)
    inverseQc = fill(NaN, nData, nCF)
    uncertaintyQc = fill(NaN, nData, nCF)
    energyRatioBCo = fill(NaN, nData, nCF)
    energyRatioCNo = fill(NaN, nData, nCF)
    raysPlot = zeros(100, 5, nData)
    tCoda = zeros(nData, nCF)
    inversionMatrixPD = zeros(nData, nModel)
    inversionMatrixQ = zeros(nData, nModel)
    inversionMatrixQc = zeros(nData, nModel, nCF)
    rayCrossing = zeros(nData, nModel)
    rejected = falses(nData)
    firstKernel = Ref{Any}(nothing)

    # Traces are independent: each thread writes only row i of the outputs, so no
    # locking is needed; `done` is an atomic counter for the progress message
    s.kernelType in ("differential", "middle") ||
        error("kernelType must be \"differential\" or \"middle\", not \"$(s.kernelType)\"")
    done = Threads.Atomic{Int}(0)
    Threads.@threads for i in 1:nData
        hdr = g.sacHeader[i]
        t, data, _ = read_sac(files[i])
        srate = 1 / hdr.delta
        sp = murat_envelope(cf, data, srate, s.envelopeSmoothTime)
        # Locations, pick (sample, time, velocity), origin time and peak delay

        locDeg, locM = murat_location(s.origin, hdr)
        locationDeg[i, :] .= locDeg
        cursorPick, pktime, v = murat_picking(t, g.PTime, g.STime, s.POrS,
                                              s.averageVelocityP, s.averageVelocityS, srate, hdr, files[i])
        theo, originTime = murat_origin_time(pktime, g.originTime, v, locM, hdr, i)
        cursorPD = murat_peak_delay_check(t, cursorPick, s.maximumPeakDelay, srate)
        pd = murat_peak_delay(sp, cursorPick, srate, cursorPD)

        # Ray tracing only on the first component of each recording
        doRays = recognize_components(i, s.components)
        if doRays
            Apd, AQ, luntot, rplot, rc = murat_rays(modv, g.gridPropagation, g.pvel, locM)
            inversionMatrixPD[i, :] .= Apd
            inversionMatrixQ[i, :] .= AQ
            totalLengthRay[i] = luntot
            raysPlot[:, :, i] .= rplot
            rayCrossing[i, :] .= rc
        end

        # Coda window; with lapseTimeMethod "Peak" tC is per frequency, and the first
        # value gives the window samples and the lapse time used for Qc
        tC, cStart, cEnd = murat_coda_check(originTime, pktime, srate, s.startLapseTime,
                                            s.codaWindow, t, pd, s.lapseTimeMethod)
        cStart, cEnd, tC1 = first(cStart), first(cEnd), first(tC)
        locationM[i, :] .= locM
        theoreticalTime[i] = theo
        tCoda[i, :] .= tC

        # Reject the trace if the coda window is cut by the end of the trace or the
        # pick-origin travel time is outside [mintravel, maxtravel] (s)
        if (cEnd - cStart) < s.codaWindow * srate - 2 ||
           (pktime - originTime) > s.maxtravel || (pktime - originTime) < s.mintravel
            rejected[i] = true
        else
            Qc, unc = murat_qc(cf, sp, cStart, cEnd, tC1, srate, s.QcMeasurement)
            if doRays
                kern(T) = murat_kernels(T, locM[1:3], locM[4:6], modvQc, s.averageVelocityS,
                                        s.kernelTreshold, s.albedo, s.iExtinctionLength,
                                        s.lapseTimeMethod)
                if s.kernelType == "differential"
                    # MuRAT 4.0: kernels at the start and at the end of the coda window
                    KS, rS = kern(tC)
                    KE, rE = kern(tC .+ s.codaWindow)
                    for j in 1:nCF
                        inversionMatrixQc[i, :, j] .=
                            murat_coda_matrix_differential(modvQc, view(KS, :, j), rS,
                                                           view(KE, :, j), rE)
                    end
                else
                    # MuRAT 3: one kernel, at the middle of the coda window
                    K, r = kern(tC .+ s.codaWindow / 2)
                    for j in 1:nCF
                        inversionMatrixQc[i, :, j] .= murat_coda_matrix(modvQc, view(K, :, j), r)
                    end
                end
                # Qc matrix rows of the first trace, for the kernel figure
                i == 1 && (firstKernel[] = (A = inversionMatrixQc[i, :, :], type = s.kernelType))
            end
            # Direct-to-coda and coda-to-noise energy ratios (coda-normalization method)
            ePS, eCN = murat_body(s.bodyWindow, s.startNoise, srate, sp, cursorPick, cStart, cEnd)
            peakDelay[i, :] .= pd
            inverseQc[i, :] .= Qc
            uncertaintyQc[i, :] .= unc
            energyRatioBCo[i, :] .= ePS
            energyRatioCNo[i, :] .= eCN
        end
        # Progress every 100 traces (counts finished traces, not the index i)
        k = Threads.atomic_add!(done, 1) + 1
        mod(k, 100) == 0 && println("Waveform number is $k")
    end

    # Store per-trace results; rayCrossing becomes the total hit count per cell
    m.rays.locationsDeg = locationDeg
    m.rays.locationsM = locationM
    m.rays.theoreticalTime = theoreticalTime
    m.rays.totalLengthRay = totalLengthRay
    m.rays.raysPlot = raysPlot
    m.rays.rayCrossing = vec(sum(rayCrossing, dims = 1))
    m.PD.peakDelay = peakDelay
    m.PD.inversionMatrixPeakDelay = inversionMatrixPD
    m.Q.inversionMatrixQ = inversionMatrixQ
    m.Q.energyRatioBodyCoda = energyRatioBCo
    m.Q.energyRatioCodaNoise = energyRatioCNo
    m.Qc.inverseQc = inverseQc
    m.Qc.uncertaintyQc = uncertaintyQc
    m.Qc.inversionMatrixQc = inversionMatrixQc
    m.Qc.tCoda = tCoda
    m.Qc.firstKernel = firstKernel[]
    m.grid.modvQc = modvQc

    murat_selection(m)
    @printf("Ratio of recordings removed due to Qc or PD inputs: %.2f%%\n", 100 * sum(rejected) / nData)
    return m
end
