# WGS84 geodesics, replacing MATLAB Mapping Toolbox `distance` and `reckon`
# with Vincenty's formulae (sub-millimetre accuracy for MuRAT distances).

"WGS84 semi-major (equatorial) axis, in m."
const WGS84_A = 6378137.0
"WGS84 flattening (dimensionless)."
const WGS84_F = 1 / 298.257223563
"WGS84 semi-minor (polar) axis, in m, derived from `WGS84_A` and `WGS84_F`."
const WGS84_B = WGS84_A * (1 - WGS84_F)

"""
    geodistance(lat1, lon1, lat2, lon2) -> (distance_m, azimuth_deg)

Distance on the WGS84 ellipsoid and forward azimuth (degrees clockwise from
north, in [0, 360)), as MATLAB `distance(lat1,lon1,lat2,lon2,wgs84Ellipsoid("m"))`.
Used to place events, stations and grid nodes in metres from the grid origin.
MATLAB: Mapping Toolbox `distance` (no MuRAT counterpart); solved here with
Vincenty's inverse formula.

# Arguments
- `lat1`, `lon1`: first point (degrees).
- `lat2`, `lon2`: second point (degrees).

# Returns
- `distance_m`: geodesic distance (m); `azimuth_deg`: azimuth at point 1 (degrees).
  Coincident points return `(0.0, 0.0)`.
"""
function geodistance(lat1::Real, lon1::Real, lat2::Real, lon2::Real)
    (lat1 == lat2 && lon1 == lon2) && return (0.0, 0.0)
    a, b, f = WGS84_A, WGS84_B, WGS84_F
    # Longitude difference and reduced latitudes U (latitudes on the auxiliary sphere)
    L = deg2rad(lon2 - lon1)
    U1 = atan((1 - f) * tan(deg2rad(lat1)))
    U2 = atan((1 - f) * tan(deg2rad(lat2)))
    sinU1, cosU1 = sincos(U1)
    sinU2, cosU2 = sincos(U2)
    # Iterate on λ (longitude difference on the auxiliary sphere) until it converges;
    # σ is the angular distance on the sphere, α the azimuth of the geodesic at the equator
    λ = L
    sinσ = cosσ = σ = cos2α = cos2σm = 0.0
    for _ in 1:200
        sinλ, cosλ = sincos(λ)
        sinσ = sqrt((cosU2 * sinλ)^2 + (cosU1 * sinU2 - sinU1 * cosU2 * cosλ)^2)
        sinσ == 0 && return (0.0, 0.0)
        cosσ = sinU1 * sinU2 + cosU1 * cosU2 * cosλ
        σ = atan(sinσ, cosσ)
        sinα = cosU1 * cosU2 * sinλ / sinσ
        cos2α = 1 - sinα^2
        # cos2α == 0: geodesic along the equator
        cos2σm = cos2α != 0 ? cosσ - 2 * sinU1 * sinU2 / cos2α : 0.0
        C = f / 16 * cos2α * (4 + f * (4 - 3 * cos2α))
        λp = λ
        λ = L + (1 - C) * f * sinα * (σ + C * sinσ * (cos2σm + C * cosσ * (-1 + 2 * cos2σm^2)))
        abs(λ - λp) < 1e-13 && break
    end
    # Series coefficients A, B and correction Δσ give the ellipsoidal distance s (m)
    u2 = cos2α * (a^2 - b^2) / b^2
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    Δσ = B * sinσ * (cos2σm + B / 4 * (cosσ * (-1 + 2 * cos2σm^2) -
                                         B / 6 * cos2σm * (-3 + 4 * sinσ^2) * (-3 + 4 * cos2σm^2)))
    s = b * A * (σ - Δσ)
    # Forward azimuth at point 1, mapped to [0, 360) degrees
    sinλ, cosλ = sincos(λ)
    α1 = atan(cosU2 * sinλ, cosU1 * sinU2 - sinU1 * cosU2 * cosλ)
    return (s, mod(rad2deg(α1), 360.0))
end

"""
    georeckon(lat1, lon1, distance_m, azimuth_deg) -> (lat2, lon2)

Point at a given distance and azimuth on the WGS84 ellipsoid, as MATLAB
`reckon(lat1,lon1,d,az,wgs84Ellipsoid("m"))`. Used to convert grid points
given in metres from the origin back to latitude/longitude.
MATLAB: Mapping Toolbox `reckon` (no MuRAT counterpart); solved here with
Vincenty's direct formula.

# Arguments
- `lat1`, `lon1`: starting point (degrees).
- `distance_m`: geodesic distance (m).
- `azimuth_deg`: azimuth at the starting point (degrees clockwise from north).

# Returns
- `lat2`, `lon2`: end point (degrees). The longitude is not wrapped to [-180, 180).
"""
function georeckon(lat1::Real, lon1::Real, s::Real, az::Real)
    a, b, f = WGS84_A, WGS84_B, WGS84_F
    α1 = deg2rad(az)
    sinα1, cosα1 = sincos(α1)
    # Reduced latitude U1 and angular distance σ1 from the equator to point 1
    tanU1 = (1 - f) * tan(deg2rad(lat1))
    cosU1 = 1 / sqrt(1 + tanU1^2)
    sinU1 = tanU1 * cosU1
    σ1 = atan(tanU1, cosα1)
    sinα = cosU1 * sinα1
    cos2α = 1 - sinα^2
    u2 = cos2α * (a^2 - b^2) / b^2
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    # Iterate on the angular distance σ on the auxiliary sphere until it converges
    σ = s / (b * A)
    cos2σm = sinσ = cosσ = 0.0
    for _ in 1:200
        cos2σm = cos(2σ1 + σ)
        sinσ, cosσ = sincos(σ)
        Δσ = B * sinσ * (cos2σm + B / 4 * (cosσ * (-1 + 2 * cos2σm^2) -
                                             B / 6 * cos2σm * (-3 + 4 * sinσ^2) * (-3 + 4 * cos2σm^2)))
        σp = σ
        σ = s / (b * A) + Δσ
        abs(σ - σp) < 1e-13 && break
    end
    # End-point latitude and longitude difference L (from λ on the auxiliary sphere)
    sinσ, cosσ = sincos(σ)
    cos2σm = cos(2σ1 + σ)
    tmp = sinU1 * sinσ - cosU1 * cosσ * cosα1
    lat2 = atan(sinU1 * cosσ + cosU1 * sinσ * cosα1, (1 - f) * sqrt(sinα^2 + tmp^2))
    λ = atan(sinσ * sinα1, cosU1 * cosσ - sinU1 * sinσ * cosα1)
    C = f / 16 * cos2α * (4 + f * (4 - 3 * cos2α))
    L = λ - (1 - C) * f * sinα * (σ + C * sinσ * (cos2σm + C * cosσ * (-1 + 2 * cos2σm^2)))
    return (rad2deg(lat2), lon1 + rad2deg(L))
end

"""
    local_xy(origin_lat, origin_lon, lat, lon) -> (x_m, y_m)

Planar offsets (m) of a point from `origin` along x (east) and y (north),
computed from the WGS84 distance and azimuth (`d sin(az)`, `d cos(az)`), as in
`Murat_location` and `Murat_checks`. All inputs in degrees.
MATLAB: no single counterpart (inline code in `Murat_location.m`).
"""
function local_xy(origin_lat, origin_lon, lat, lon)
    d, az = geodistance(origin_lat, origin_lon, lat, lon)
    return d * sin(deg2rad(az)), d * cos(deg2rad(az))
end

"""
    km2deg(km) -> degrees

MATLAB `km2deg` (spherical Earth of radius 6371 km): converts a distance in km
to degrees of arc.
"""
km2deg(km) = km / (6371.0 * π / 180)
