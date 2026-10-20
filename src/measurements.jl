# Envelopes and per-trace measurements: picking, origin time, peak delay,
# coda window, Qc and body/coda/noise energy ratios.

"""
    murat_envelope(cf, data, srate, envelopeSmoothTime) -> sp

`Murat_envelope`: for each central frequency, 4th-order Butterworth
band-pass (cf ± cf/3, zero phase), then the squared RMS envelope over a sliding
window of `envelopeSmoothTime` seconds. Returns a matrix (samples × frequencies).
These energy envelopes are the input of all the per-trace measurements.
MATLAB: `Murat_envelope.m` of MuRAT 4.0, `envelope(fsisma, winSamples, 'rms').^2`
(the trace is read by the caller here, not inside). MuRAT 3 used the squared
Hilbert envelope smoothed over 6 periods instead.

# Arguments
- `cf`: central frequencies (Hz).
- `data`: seismogram samples (one trace).
- `srate`: sampling rate (Hz).
- `envelopeSmoothTime`: length of the RMS window (s); the window has
  `max(1, round(envelopeSmoothTime * srate))` samples.

# Returns
- `sp`: squared RMS envelopes, length(data) × length(cf). A column stays
  zero if its band-pass is invalid (a warning is printed).
"""
function murat_envelope(cf::AbstractVector, data::AbstractVector, srate::Real,
                        envelopeSmoothTime::Real)
    (srate == -12345 || !isfinite(srate) || srate <= 0) && error("Waveform has no valid sampling rate!")
    # Remove the linear trend and taper 5% of the trace (Tukey window)
    sisma = detrend_linear(data)
    lsis = length(sisma)
    tsisma = tukeywin(lsis, 0.05) .* sisma
    sp = zeros(lsis, length(cf))
    # RMS window in samples (e.g. 0.25 s at 100 Hz = 25 samples)
    winSamples = max(1, round(Int, envelopeSmoothTime * srate, RoundNearestTiesAway))
    for (i, f) in enumerate(cf)
        # Band corners normalized to the Nyquist frequency, clipped inside (0, 1)
        Wn = [f - f / 3, f + f / 3] ./ srate .* 2
        Wn[1] <= 0 && (Wn[1] = 0.001)
        Wn[2] >= 1 && (Wn[2] = 0.999)
        if Wn[2] <= Wn[1]
            @warn "Skipping frequency $f Hz: invalid bandpass."
            continue
        end
        # Zero-phase filtering in second-order sections (MATLAB zp2sos + filtfilt);
        # energy = squared RMS envelope
        filt = digitalfilter(Bandpass(Wn[1], Wn[2]), Butterworth(4))
        fsisma = filtfilt(convert(SecondOrderSections, filt), tsisma)
        sp[:, i] = rms_envelope(fsisma, winSamples) .^ 2
    end
    return sp
end

"""
    murat_location(origin, hdr) -> (locDeg, locM)

`Murat_location`: event and station coordinates in degrees
`[evla evlo evdp(m) stla stlo stel]` and in meters from the grid origin
`[x y z]` of event and station.
The metric coordinates are used for ray tracing and kernels; the degrees for
plots and declustering.
MATLAB: `Murat_location.m`.

# Arguments
- `origin`: grid origin `[lat, lon, altitude (m)]` (degrees, degrees, m).
- `hdr`: SAC header of the trace (`evdp` in km, `stel` in m).

# Returns
- `locDeg`: 6-vector `[evla, evlo, evdp, stla, stlo, stel]` (degrees, depth/elevation in m).
- `locM`: 6-vector `[xe, ye, ze, xs, ys, zs]` (m; x east, y north, z up, so the
  event z is minus its depth).
"""
function murat_location(origin, hdr::SacHeader)
    # Event depth from km to m; SAC_UNDEF (-12345) marks a missing header field
    even = (hdr.evla, hdr.evlo, hdr.evdp * 1000)
    stati = (hdr.stla, hdr.stlo, hdr.stel)
    (any(==(SAC_UNDEF), even) || any(==(SAC_UNDEF), stati)) &&
        error("Waveform has missing event or station location fields!")
    locDeg = [even..., stati...]
    xe, ye = local_xy(origin[1], origin[2], even[1], even[2])
    xs, ys = local_xy(origin[1], origin[2], stati[1], stati[2])
    locM = [xe, ye, -even[3], xs, ys, stati[3]]
    return locDeg, locM
end

"""
    murat_picking(t, PTime, STime, PorS, vP, vS, srate, hdr, name)
        -> (cursorPick, pktime, v)

`Murat_picking`: pick time, its sample and the average velocity of the phase.
The pick is read from the SAC header field named by `PTime` (P, `PorS == 2`) or
`STime` (S, `PorS == 3`); it must lie inside the trace.
MATLAB: `Murat_picking.m` (MATLAB evaluates the header expression with `eval`).

# Arguments
- `t`: time vector of the trace (s, SAC reference time).
- `PTime`, `STime`: header field names of the P and S picks (e.g. `"a"`, `"t0"`).
- `PorS`: 2 for P, 3 for S.
- `vP`, `vS`: average P and S velocities (km/s).
- `srate`: sampling rate (Hz); `hdr`: SAC header; `name`: file name (for errors).

# Returns
- `cursorPick`: sample offset of the pick from the trace start, `floor((pk - t[1]) srate)`
  (used directly as a 1-based index, as in MATLAB).
- `pktime`: pick time (s); `v`: average velocity of the picked phase (km/s).
"""
function murat_picking(t, PTime, STime, PorS, vP, vS, srate, hdr, name)
    if PorS == 2
        pk, v = sac_time(hdr, PTime), vP
    elseif PorS == 3
        pk, v = sac_time(hdr, STime), vS
    else
        error("POrS must be 2 (P) or 3 (S)")
    end
    pk < t[1] && error("Picking is before the start of the waveform $name")
    pk > t[end] && error("Picking is after the end of the waveform $name")
    return floor(Int, (pk - t[1]) * srate), pk, v
end

"""
    murat_origin_time(pk, originTime, v, locM, hdr, i) -> (theo, origin)

`Murat_originTime`: origin time from the header or from the straight-ray travel time.
If the header origin time is undefined, the travel time is the straight
event-station distance divided by the average velocity, and the origin time
is the pick minus it.
MATLAB: `Murat_originTime.m`.

# Arguments
- `pk`: pick time (s); `originTime`: header field name of the origin time.
- `v`: average velocity of the phase (km/s); `locM`: `[xe ye ze xs ys zs]` (m).
- `hdr`: SAC header; `i`: trace number (for errors).

# Returns
- `theo`: travel time from origin to pick (s); `origin`: origin time (s).
"""
function murat_origin_time(pk, originTime, v, locM, hdr, i)
    o = sac_time(hdr, originTime)
    if o == SAC_UNDEF
        dist = sqrt((locM[4] - locM[1])^2 + (locM[5] - locM[2])^2 + (locM[6] - locM[3])^2)
        # Distance in m, velocity in km/s -> travel time in s
        theo = dist / v / 1000
        return theo, pk - theo
    end
    theo = pk - o
    theo < 0 && error("Picking set before the origin time, recording $i")
    return theo, o
end

"""
    murat_peak_delay_check(t, cursorPick, maxtpde, srate) -> cursorPeakDelay

`Murat_peakDelayCheck`: last sample of the peak-delay window.
The window spans `maxtpde` seconds after the pick and is cut at the end of the
trace (with a warning).
MATLAB: `Murat_peakDelayCheck.m`.

# Arguments
- `t`: time vector (s); `cursorPick`: pick sample; `maxtpde`: maximum peak delay (s);
  `srate`: sampling rate (Hz).
"""
function murat_peak_delay_check(t, cursorPick, maxtpde, srate)
    cmax = floor(Int, cursorPick + maxtpde * srate)
    c = min(cmax, length(t))
    c < cmax && @warn "Length of seismogram shorter than peak-delay window!"
    return c
end

"""
    murat_peak_delay(sp, cursorPick, srate, cursorPeakDelay) -> peakDelay

`Murat_peakDelay`: time of the envelope maximum after the pick, per frequency.
The peak delay (time from the direct-wave onset to the envelope maximum)
measures scattering along the path.
MATLAB: `Murat_peakDelay.m`.

# Arguments
- `sp`: envelopes (samples × frequencies); `cursorPick`, `cursorPeakDelay`: first
  and last samples of the search window; `srate`: sampling rate (Hz).

# Returns
- `peakDelay`: vector (one per frequency, s). The index of the maximum counts from
  1 at the pick sample, as in MATLAB, so the value is one sample longer than the
  offset from the pick.
"""
function murat_peak_delay(sp, cursorPick, srate, cursorPeakDelay)
    w = @view sp[cursorPick:cursorPeakDelay, :]
    return [argmax(view(w, :, j)) for j in axes(w, 2)] ./ srate
end

"""
    murat_coda_check(originTime, pktime, srate, tCm, tWm, t, peakDelay, method)
        -> (tCoda, cStart, cEnd)

`Murat_codaCheck`: lapse time of the coda start and its first and last
samples on the trace.
The coda starts after the envelope peak (`"Peak"`), at a fixed lapse time
(`"Constant"`) or at a multiple of the travel time (`"Travel"`); the window
lasts `tWm` seconds and is cut at the end of the trace.
MATLAB: `Murat_codaCheck.m`.

# Arguments
- `originTime`, `pktime`: origin and pick times (s).
- `srate`: sampling rate (Hz); `t`: time vector (s).
- `tCm`: coda start lapse time (s, `"Constant"`) or travel-time factor (`"Travel"`).
- `tWm`: coda window length (s).
- `peakDelay`: peak delays (s, one per frequency), used by `"Peak"`.
- `method`: `"Peak"`, `"Constant"` or `"Travel"` (input `lapseTimeMethod`).

# Returns
- `tCoda`: lapse time of the coda start from the origin (s); a vector (one per
  frequency) for `"Peak"`, a scalar otherwise.
- `cStart`, `cEnd`: first and last coda samples (same shape as `tCoda`).
"""
function murat_coda_check(originTime, pktime, srate, tCm, tWm, t, peakDelay, method)
    pk_rel = pktime - originTime
    tCoda = if method == "Peak"
        pk_rel .+ peakDelay
    elseif method == "Constant"
        tCm
    elseif method == "Travel"
        originTime + tCm * pk_rel
    else
        error("Unknown peak-delay method")
    end
    # Lapse time -> absolute time -> sample from the trace start (the -1 is from MATLAB)
    cStart = floor.(Int, (tCoda .+ originTime .- t[1]) .* srate .- 1)
    cEnd = min.(length(t), cStart .+ floor(Int, tWm * srate - 1))
    return tCoda, cStart, cEnd
end

"""
Trial inverse-Qc values of the `"NonLinear"` grid search in `murat_qc`
(MATLAB `10^-5:10^-5:10^-1`, built with the MATLAB colon rules).
"""
const QVALUES_NONLINEAR = matlab_colon(1e-5, 1e-5, 1e-1)

"""
    murat_qc(cf, sp, cStart, cEnd, tCoda, srate, method) -> (Qc, unc)

`Murat_Qc`: inverse coda attenuation and its uncertainty per frequency, by
linear regression (`"Linearized"`) or grid search on 1-s averages
(`"NonLinear"`).
Both fit the single-scattering coda model E(t) ∝ t^-1.5 exp(-2π f t Qc^-1)
to the coda envelope.
MATLAB: `Murat_Qc.m` (`polyfit` / `corrcoef` replaced by a least-squares solve
and `cor`).

# Arguments
- `cf`: central frequencies (Hz); `sp`: envelopes (samples × frequencies).
- `cStart`, `cEnd`: first and last coda samples; `tCoda`: lapse time of the coda
  start (s); `srate`: sampling rate (Hz, rounded to an integer as in MATLAB).
- `method`: `"Linearized"` or `"NonLinear"` (input `QcMeasurement`).

# Returns
- `Qc`: inverse coda Q per frequency (0 where the fit fails).
- `unc`: uncertainty per frequency: `1/|R|` (R = correlation coefficient of the
  linear fit) for `"Linearized"`, a misfit-based value for `"NonLinear"`; 0 where
  `Qc` is 0.
"""
function murat_qc(cf, sp, cStart, cEnd, tCoda, srate, method)
    lcf = length(cf)
    Qc = zeros(lcf)
    unc = zeros(lcf)
    srate = round(Int, srate, RoundNearestTiesAway)
    si = 1 / srate
    for i in 1:lcf
        env = sp[cStart:cEnd, i]
        lenv = length(env)
        cf_i = cf[i]
        if method == "Linearized"
            # Exactly lenv samples: the float colon can drop one (tCoda+si:si:tCoda+lenv*si)
            lapseT = tCoda .+ (1:lenv) .* si
            # Drop 5% of the samples at each end; fit log(E t^1.5)/(2π f) = -Qc^-1 t + c,
            # so the slope is minus the inverse Qc
            edgeC = floor(Int, 0.05 * lenv)
            lt = lapseT[edgeC:end-edgeC]
            e = env[edgeC:end-edgeC]
            y = log.(e .* lt .^ 1.5) ./ (2π * cf_i)
            G = hcat(lt, ones(length(lt)))
            p = G \ y
            R = cor(lt, y)
            if p[1] < 0
                Qc[i] = -p[1]
                unc[i] = 1 / abs(R)
            end
        elseif method == "NonLinear"
            # Average the envelope over consecutive 1-s blocks (srate samples each)
            lWindow = round(Int, (cEnd - cStart) * si, RoundNearestTiesAway)
            nCols = min(lWindow, fld(length(env), srate))
            blk = nCols == 0 ? zeros(srate, 0) : reshape(env[1:nCols*srate], srate, nCols)
            dObs = vec(mean(blk, dims = 1))
            # Normalize by the last block (removes the unknown amplitude); L = block
            # centres in lapse time (s). L1 misfit for each trial inverse Qc
            dObserved = dObs[1:end-1] ./ dObs[end]
            L = matlab_colon(tCoda + 0.5, 1, tCoda + size(blk, 2) - 0.5)
            c = 2π * cf_i
            E = zeros(length(QVALUES_NONLINEAR))
            for (q, Qv) in enumerate(QVALUES_NONLINEAR)
                pred = L .^ -1.5 .* exp.(-c .* L .* Qv)
                E[q] = sum(abs.(dObserved .- pred[1:end-1] ./ pred[end]))
            end
            Emin, idx = findmin(E)
            fitQ = QVALUES_NONLINEAR[idx]
            # Best fit at the lower bound of the search counts as no measurement
            if fitQ != minimum(QVALUES_NONLINEAR)
                mL = mean(L)
                Qc[i] = fitQ
                unc[i] = abs(log(Emin / mL^1.5) / (c * mL))
            end
        else
            error("Unknown strategy to calculate Qc")
        end
    end
    # Force zero uncertainty where Qc is zero (as in MATLAB)
    unc[Qc .<= 0] .= 0
    Qc[unc .< 0] .= 0
    return Qc, unc
end

"""
    murat_body(bodyWindow, startNoise, srate, sp, cursorPick, cStart, cEnd)
        -> (ratioBodyCoda, ratioCodaNoise)

`Murat_body`: direct-to-coda and coda-to-noise energy ratios per frequency.
Energies are time integrals of the envelope in the direct-wave window (from the
pick), in the noise window and in the coda window; the coda energy is scaled to
the length of the direct-wave window. These ratios feed the coda-normalization
(Q) inversion and its noise threshold.
MATLAB: `Murat_body.m`.

# Arguments
- `bodyWindow`: direct-wave window length (s).
- `startNoise`: start of the noise window (s from the trace start).
- `srate`: sampling rate (Hz); `sp`: envelopes (samples × frequencies).
- `cursorPick`: pick sample; `cStart`, `cEnd`: first and last coda samples.

# Returns
- `ratioBodyCoda`, `ratioCodaNoise`: vectors, one value per frequency.
"""
function murat_body(bodyWindow, startNoise, srate, sp, cursorPick, cStart, cEnd)
    # Window bounds in samples; the noise window has the direct-wave window length
    intSamples = floor(Int, bodyWindow * srate)
    n0 = floor(Int, startNoise * srate)
    n1 = n0 + intSamples
    ps1 = floor(Int, cursorPick + intSamples - 1)
    lengthSamples = length(cStart:cEnd)
    dt = 1 / srate
    # Trapezoidal integrals; coda energy normalized by (coda length / direct length)
    PS = trapz_columns(sp[cursorPick:ps1, :], dt)
    noise = trapz_columns(sp[n0:n1, :], dt)
    coda = trapz_columns(sp[cStart:cEnd, :], dt) ./ (lengthSamples / intSamples)
    return PS ./ coda, coda ./ noise
end
