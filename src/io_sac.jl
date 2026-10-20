# SAC file reading (replaces MatSAC fget_sac/sac/sachdr) and data checks.

# Value of undefined SAC header fields
const SAC_UNDEF = -12345.0

"""
    SacHeader

Header values of a SAC file that MuRAT uses (same names as MatSAC).
MATLAB: the `SAChdr` structure returned by MatSAC `fget_sac` (e.g. `SAChdr.times.a`,
`SAChdr.station.stla`, `SAChdr.event.evla`), flattened into one immutable struct.

# Fields
- `delta`: sampling interval (s); `b`, `e`: first and last sample times (s);
  `o`: origin time; `a`: first arrival (P) pick; `t`: picks `t0` ... `t9` (s).
- `stla`, `stlo`, `stel`: station latitude, longitude (degrees), elevation (m).
- `evla`, `evlo`, `evdp`: event latitude, longitude (degrees), depth (as in the file).
- `mag`: magnitude; `npts`: number of samples; `kstnm`: station name;
  `kcmpnm`: component name. Undefined values are `SAC_UNDEF` (-12345).
"""
struct SacHeader
    delta::Float64
    b::Float64
    e::Float64
    o::Float64
    a::Float64
    t::NTuple{10,Float64}     # t0 ... t9
    stla::Float64
    stlo::Float64
    stel::Float64
    evla::Float64
    evlo::Float64
    evdp::Float64
    mag::Float64
    npts::Int
    kstnm::String
    kcmpnm::String
end

"""
    sac_time(hdr, name)

Header time selected by name, as in the MuRAT input (`"o"`, `"a"`, `"t0"`
... `"t9"`), also accepting the MATLAB form `"SAChdr.times.t0"`.
`"b"` and `"e"` are accepted too. Returns the time in s (`SAC_UNDEF` if undefined).
MATLAB: `eval(originTime)` / `getFieldValue` in `Murat_testData.m` and the
`SAChdr.times.<name>` reads in the data loop.
"""
function sac_time(hdr::SacHeader, name::AbstractString)
    n = replace(name, "SAChdr.times." => "")
    n == "o" && return hdr.o
    n == "a" && return hdr.a
    n == "b" && return hdr.b
    n == "e" && return hdr.e
    # "tK" -> hdr.t[K + 1] (t0 is stored first)
    if startswith(n, "t") && length(n) == 2 && isdigit(n[2])
        return hdr.t[parse(Int, n[2:2])+1]
    end
    error("Unknown SAC time field '$name'")
end

"""
    read_sac(filename) -> (t, data, hdr)

Reads a binary SAC file (either byte order). `t` is the time vector
`b + (0:npts-1)*delta`, `data` the samples as Float64.
MATLAB: MatSAC `fget_sac` (called through `Murat_test.m`). Only evenly sampled,
single-component time series are supported (no second data block is read).

# Arguments
- `filename`: path of the SAC file.

# Returns
- `t`: sample times (s, relative to the SAC reference time), length `npts`.
- `data`: samples (Float64), length `npts`.
- `hdr`: `SacHeader`.
"""
function read_sac(filename::AbstractString)
    open(filename, "r") do io
        # SAC header: 70 floats, 40 integers, 192 characters (632 bytes), then data
        fl = read!(io, Vector{Float32}(undef, 70))
        it = read!(io, Vector{Int32}(undef, 40))
        ch = read!(io, Vector{UInt8}(undef, 192))
        # Byte order: the header version (integer word 7) is small only if read in the
        # right order
        swap = !(1 <= it[7] <= 20)          # nvhdr is 6 or 7
        if swap
            fl = bswap.(fl)
            it = bswap.(it)
        end
        npts = Int(it[10])
        data = read!(io, Vector{Float32}(undef, npts))
        swap && (data = bswap.(data))
        # Closures: float header word i (1-based) and trimmed character field r
        f(i) = Float64(fl[i])
        str(r) = strip(String(ch[r]), ['\0', ' '])
        # 1-based float words: delta 1, b 6, e 7, o 8, a 9, t0-t9 11-20, stla 32,
        # stlo 33, stel 34, evla 36, evlo 37, evdp 39, mag 40; characters: kstnm
        # 1-8, kcmpnm 161-168
        hdr = SacHeader(f(1), f(6), f(7), f(8), f(9), ntuple(k -> f(10 + k), 10),
                        f(32), f(33), f(34), f(36), f(37), f(39), f(40), npts,
                        str(1:8), str(161:168))
        t = hdr.b .+ (0:npts-1) .* hdr.delta
        return collect(t), Float64.(data), hdr
    end
end

"""
    list_sac(folder) -> Vector{String}

Lists the SAC files of a folder (any extension case, hidden files skipped).
Returns full paths sorted by file name.
MATLAB: `createsList([dataDirectory '/*.sac'])` in `Murat_checks.m`. Difference:
MATLAB matches only `.sac` (case-sensitive on Linux), here `.SAC` is also accepted.
"""
function list_sac(folder::AbstractString)
    names = filter(n -> !startswith(n, ".") && lowercase(splitext(n)[2]) == ".sac", readdir(folder))
    return [joinpath(folder, n) for n in sort(names)]
end

"""
    test_data(files, originTime, PTime, STime)

`Murat_testData`: reads all headers, checks the mandatory fields and
returns the headers, a table of the main values and a flag for missing
origin (1) or S (2) times.
MATLAB: `Murat_testData.m`. Errors if a P time or an event/station coordinate is
missing; missing origin or S times are only flagged.

# Arguments
- `files`: SAC file paths.
- `originTime`, `PTime`, `STime`: header field names of the three times (`"o"`, ...).

# Returns
- `headers`: `Vector{SacHeader}`, one per file.
- `table`: vector of NamedTuples (`Names`, `Origin`, `P`, `S`, `EvLat`, `EvLon`,
  `EvDepth`, `StLat`, `StLon`, `StElev`), undefined times as `missing`.
- `flag`: 0 = all present, 1 = origin missing, 2 = S missing, 3 = both missing
  (in at least one file).
"""
function test_data(files, originTime, PTime, STime)
    headers = Vector{SacHeader}(undef, length(files))
    table = Vector{NamedTuple}(undef, length(files))
    origMissing = false
    sMissing = false
    for (i, f) in enumerate(files)
        _, _, h = read_sac(f)
        headers[i] = h
        o = sac_time(h, originTime)
        p = sac_time(h, PTime)
        s = sac_time(h, STime)
        o == SAC_UNDEF && (origMissing = true)
        p == SAC_UNDEF && error("There is a missing P value, check $f!")
        s == SAC_UNDEF && (sMissing = true)
        for (v, what) in ((h.evla, "event"), (h.evlo, "event"), (h.evdp, "event"),
                          (h.stla, "station"), (h.stlo, "station"), (h.stel, "station"))
            v == SAC_UNDEF && error("There is a missing $what coordinate, check $f!")
        end
        table[i] = (Names = f, Origin = o == SAC_UNDEF ? missing : o, P = p,
                    S = s == SAC_UNDEF ? missing : s, EvLat = h.evla, EvLon = h.evlo,
                    EvDepth = h.evdp, StLat = h.stla, StLon = h.stlo, StElev = h.stel)
    end
    # Bool arithmetic: origin -> 1, S -> 2, both -> 3 (as MATLAB origMissing + 2*sMissing)
    return headers, table, origMissing + 2sMissing
end
