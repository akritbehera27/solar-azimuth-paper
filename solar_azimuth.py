#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
solar_azimuth.py
================
"Where does the Sun rise?" -- a solar-azimuth engine for the 22 km region
centred on  23.169754 N, 82.360300 E  (Khadganva, Manendragarh-Chirimiri-
Bharatpur district, Chhattisgarh, India).

The file contains THREE independent models of the same sky, on purpose:

  Level 2  "noaa_full"   : NOAA's operational JavaScript algorithm
                           (gml.noaa.gov/grad/solcalc/main.js), transcribed
                           line-by-line.  Used as the accuracy benchmark.
  Level 1  "noaa_pdf"    : the closed-form NOAA sheet
                           (gml.noaa.gov/grad/solcalc/solareqns.PDF) --
                           fractional year, equation of time, declination,
                           hour angle, zenith, azimuth.  The "verified source".
  Level 0  "simple"      : the Class-12 single-harmonic model
                           (delta = 23.44 sin(...), EoT = 9.87 sin2B - ...),
                           plus the refraction-corrected sunrise formula.

All azimuths are measured CLOCKWISE FROM TRUE NORTH:

    0 deg = N, 90 deg = E, 180 deg = S, 270 deg = W.

Run:
    python3 solar_azimuth.py verify     # cross-check the models (needs node)
    python3 solar_azimuth.py tables     # write the CSV tables used in the paper
    python3 solar_azimuth.py demo       # print the worked example, step by step
"""
from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
from datetime import date, timedelta

# ----------------------------------------------------------------------
# geometry helpers (degrees in, degrees out -- friendlier for students)
# ----------------------------------------------------------------------
def sind(x): return math.sin(math.radians(x))
def cosd(x): return math.cos(math.radians(x))
def tand(x): return math.tan(math.radians(x))
def asind(x): return math.degrees(math.asin(max(-1.0, min(1.0, x))))
def acosd(x): return math.degrees(math.acos(max(-1.0, min(1.0, x))))
def atan2d(y, x): return math.degrees(math.atan2(y, x))
def wrap360(a): return a % 360.0
def wrap180(a):
    a = a % 360.0
    return a - 360.0 if a > 180.0 else a


# ----------------------------------------------------------------------
# the site and its 22 km region
# ----------------------------------------------------------------------
SITE = {
    "name": "Khadganva (Manendragarh-Chirimiri-Bharatpur dist.), Chhattisgarh",
    "lat": 23.169753993380713,      # degrees North
    "lon": 82.36029982215409,       # degrees East
    "tz": 5.5,                      # IST = UTC + 5:30
    "tz_name": "IST (UTC+05:30)",
    "lon_std": 82.5,                # standard meridian of IST
}
REGION_RADIUS_KM = 22.0
KM_PER_DEG_LAT = 111.32             # mean; true value 110.57 at the equator -> 111.69 at the poles


def region_bounds(lat, lon, radius_km=REGION_RADIUS_KM):
    """Lat/lon half-widths of a circle of `radius_km` around (lat, lon)."""
    dlat = radius_km / KM_PER_DEG_LAT
    dlon = radius_km / (KM_PER_DEG_LAT * cosd(lat))
    return dlat, dlon


# ----------------------------------------------------------------------
# calendar
# ----------------------------------------------------------------------
def is_leap(y): return (y % 4 == 0 and y % 100 != 0) or y % 400 == 0
def days_in_year(y): return 366 if is_leap(y) else 365


def day_of_year(y, m, d):
    """Days elapsed since 1 January (1 Jan = 1)."""
    days_before = [0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    n = days_before[m] + d
    if m > 2 and is_leap(y):
        n += 1
    return n


def jd_from_ymd(y, m, d):
    """Julian Day at 00:00 UT (NOAA getJD)."""
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    return math.floor(365.25 * (y + 4716)) + math.floor(30.6001 * (m + 1)) + d + B - 1524.5


def julian_centuries(jd):
    return (jd - 2451545.0) / 36525.0


# ======================================================================
# LEVEL 2 -- NOAA's full published JS algorithm, transcribed exactly
# ======================================================================
def _L0(t):  return wrap360(280.46646 + t * (36000.76983 + t * 0.0003032))
def _M(t):   return 357.52911 + t * (35999.05029 - 0.0001537 * t)
def _ecc(t): return 0.016708634 - t * (0.000042037 + 0.0000001267 * t)


def _sun_eq_center(t):
    m = _M(t)
    return (sind(m) * (1.914602 - t * (0.004817 + 0.000014 * t))
            + sind(2 * m) * (0.019993 - 0.000101 * t)
            + sind(3 * m) * 0.000289)


def _apparent_long(t):
    o = _L0(t) + _sun_eq_center(t)
    omega = 125.04 - 1934.136 * t
    return o - 0.00569 - 0.00478 * sind(omega)


def _obliquity(t):
    secs = 21.448 - t * (46.8150 + t * (0.00059 - t * 0.001813))
    e0 = 23.0 + (26.0 + secs / 60.0) / 60.0
    return e0 + 0.00256 * cosd(125.04 - 1934.136 * t)


def noaa_full_decl(t):
    """Solar declination, degrees (NOAA main.js)."""
    return asind(sind(_obliquity(t)) * sind(_apparent_long(t)))


def noaa_full_radius(t):
    """Sun-Earth distance in AU (NOAA main.js)."""
    v = _M(t) + _sun_eq_center(t)
    e = _ecc(t)
    return (1.000001018 * (1 - e * e)) / (1 + e * cosd(v))


def noaa_full_eqtime(t):
    """Equation of time, minutes (NOAA main.js)."""
    eps, l0, e, m = _obliquity(t), _L0(t), _ecc(t), _M(t)
    y = tand(eps / 2.0) ** 2
    E = (y * sind(2 * l0)
         - 2.0 * e * sind(m)
         + 4.0 * e * y * sind(m) * cosd(2 * l0)
         - 0.5 * y * y * sind(4 * l0)
         - 1.25 * e * e * sind(2 * m))
    return 4.0 * math.degrees(E)


def noaa_refraction(elev_deg):
    """Atmospheric refraction correction in degrees (NOAA main.js)."""
    if elev_deg > 85.0:
        return 0.0
    te = tand(elev_deg)
    if elev_deg > 5.0:
        corr = 58.1 / te - 0.07 / te ** 3 + 0.000086 / te ** 5
    elif elev_deg > -0.575:
        corr = 1735.0 + elev_deg * (-518.2 + elev_deg * (103.4 + elev_deg * (-12.79 + elev_deg * 0.711)))
    else:
        corr = -20.774 / te
    return corr / 3600.0


def noaa_full_azel(t, minutes_local, lat, lon, tz):
    """Azimuth & elevation (deg) at a clock time, NOAA main.js calcAzEl()."""
    eqtime = noaa_full_eqtime(t)
    decl = noaa_full_decl(t)
    tst = minutes_local + eqtime + 4.0 * lon - 60.0 * tz
    while tst > 1440.0:
        tst -= 1440.0
    ha = tst / 4.0 - 180.0
    if ha < -180.0:
        ha += 360.0
    csz = sind(lat) * sind(decl) + cosd(lat) * cosd(decl) * cosd(ha)
    zenith = acosd(csz)
    denom = cosd(lat) * sind(zenith)
    if abs(denom) > 0.001:
        az_rad = (sind(lat) * cosd(zenith) - sind(decl)) / denom
        azimuth = 180.0 - acosd(max(-1.0, min(1.0, az_rad)))
        if ha > 0.0:
            azimuth = -azimuth
    else:
        azimuth = 180.0 if lat > 0.0 else 0.0
    azimuth = wrap360(azimuth)
    elevation = 90.0 - (zenith - noaa_refraction(90.0 - zenith))
    return azimuth, elevation, ha, decl, eqtime, zenith


def noaa_full_sunrise_hourangle(lat, decl, z0=90.833):
    """Hour angle at sunrise/sunset, radians (NOAA main.js calcHourAngleSunrise)."""
    arg = cosd(z0) / (cosd(lat) * cosd(decl)) - tand(lat) * tand(decl)
    return math.acos(max(-1.0, min(1.0, arg)))


def _noaa_rise_set_utc(rise, jd, lat, lon):
    t = julian_centuries(jd)
    eqtime = noaa_full_eqtime(t)
    decl = noaa_full_decl(t)
    ha = noaa_full_sunrise_hourangle(lat, decl)
    if not rise:
        ha = -ha
    delta = lon + math.degrees(ha)
    return 720.0 - 4.0 * delta - eqtime


def noaa_full_rise_set(rise, y, m, d, lat, lon, tz):
    """NOAA calcSunriseSet(): returns (local minutes, azimuth, declination)."""
    jd = jd_from_ymd(y, m, d)
    t_utc = _noaa_rise_set_utc(rise, jd, lat, lon)
    new_utc = _noaa_rise_set_utc(rise, jd + t_utc / 1440.0, lat, lon)
    t_local = new_utc + tz * 60.0
    t = julian_centuries(jd + new_utc / 1440.0)
    az, _, _, decl, _, _ = noaa_full_azel(t, t_local, lat, lon, tz)
    return t_local, az, decl


def noaa_full_solar_noon(y, m, d, lon, tz):
    jd = jd_from_ymd(y, m, d)
    t = julian_centuries(jd - lon / 360.0)
    off = 720.0 - lon * 4.0 - noaa_full_eqtime(t)
    t2 = julian_centuries(jd - 0.5 + off / 1440.0)
    noon = 720.0 - lon * 4.0 - noaa_full_eqtime(t2) + tz * 60.0
    return noon % 1440.0


# ======================================================================
# LEVEL 1 -- the NOAA closed-form sheet (solareqns.PDF)
# ======================================================================
def noaa_pdf_gamma(doy, hour):
    return (2.0 * math.pi / 365.0) * (doy - 1.0 + (hour - 12.0) / 24.0)


def noaa_pdf_eqtime(gamma):
    return 229.18 * (0.000075
                     + 0.001868 * math.cos(gamma) - 0.032077 * math.sin(gamma)
                     - 0.014615 * math.cos(2 * gamma) - 0.040849 * math.sin(2 * gamma))


def noaa_pdf_decl(gamma):
    return math.degrees(0.006918
                        - 0.399912 * math.cos(gamma) + 0.070257 * math.sin(gamma)
                        - 0.006758 * math.cos(2 * gamma) + 0.000907 * math.sin(2 * gamma)
                        - 0.002697 * math.cos(3 * gamma) + 0.00148 * math.sin(3 * gamma))


def noaa_pdf_azel(y, m, d, minutes_local, lat, lon, tz):
    """Azimuth (deg, from North, clockwise), zenith and hour angle, NOAA PDF route."""
    doy = day_of_year(y, m, d)
    hr = minutes_local / 60.0
    gamma = noaa_pdf_gamma(doy, hr)
    eqtime = noaa_pdf_eqtime(gamma)
    decl = noaa_pdf_decl(gamma)
    tst = minutes_local + eqtime + 4.0 * lon - 60.0 * tz
    ha = wrap180(tst / 4.0 - 180.0)
    csz = sind(lat) * sind(decl) + cosd(lat) * cosd(decl) * cosd(ha)
    zenith = acosd(csz)
    denom = cosd(lat) * sind(zenith)
    if abs(denom) > 0.001:
        az = 180.0 - acosd(max(-1.0, min(1.0, (sind(lat) * cosd(zenith) - sind(decl)) / denom)))
        if ha > 0.0:
            az = -az
    else:
        az = 180.0 if lat > 0 else 0.0
    return wrap360(az), zenith, ha, decl, eqtime


# ======================================================================
# LEVEL 0 -- the Class-12 model
# ======================================================================
def simple_decl(n):
    """Declination, degrees: 23.44 sin(360/365 (N-81))."""
    return 23.44 * sind(360.0 / 365.0 * (n - 81.0))


def simple_eot(n):
    """Equation of time, minutes: 9.87 sin2B - 7.53 cosB - 1.5 sinB."""
    B = 360.0 * (n - 81.0) / 364.0
    return 9.87 * sind(2 * B) - 7.53 * cosd(B) - 1.5 * sind(B)


def clock_to_solar_time(clock_hours, n, lon, lon_std=82.5):
    """IST clock time -> local apparent solar time (hours)."""
    return clock_hours + simple_eot(n) / 60.0 + (lon - lon_std) / 15.0


def solar_time_to_clock(solar_hours, n, lon, lon_std=82.5):
    return solar_hours - simple_eot(n) / 60.0 - (lon - lon_std) / 15.0


def hour_angle(solar_hours):
    return 15.0 * (solar_hours - 12.0)


def zenith_angle(lat, decl, H):
    """z from the spherical triangle: cos z = sin f sin d + cos f cos d cos H."""
    return acosd(sind(lat) * sind(decl) + cosd(lat) * cosd(decl) * cosd(H))


def azimuth_general(lat, decl, H):
    """Azimuth (deg, from North, clockwise) for hour angle H (deg).

    A = 180 + atan2( sin H , cos H sin(lat) - tan(decl) cos(lat) )
    -- one equation, no AM/PM case-work, valid for the whole day."""
    return wrap360(180.0 + atan2d(sind(H), cosd(H) * sind(lat) - tand(decl) * cosd(lat)))


def azimuth_arccos(lat, decl, H):
    """Azimuth from the arccos form of the NOAA sheet (needs the AM/PM rule)."""
    z = zenith_angle(lat, decl, H)
    az = 180.0 - acosd((sind(lat) * cosd(z) - sind(decl)) / (cosd(lat) * sind(z)))
    if H > 0.0:
        az = -az
    return wrap360(az)


def sunrise_hour_angle(lat, decl, z0=90.833):
    """H0 in degrees: cos H0 = cos z0/(cos f cos d) - tan f tan d."""
    return acosd(cosd(z0) / (cosd(lat) * cosd(decl)) - tand(lat) * tand(decl))


def sunrise_azimuth(lat, decl, z0=90.833):
    """Sunrise azimuth (deg from North, clockwise) on a day of declination `decl`.

        cos A = [ sin(d) - sin(f) cos(z0) ] / [ cos(f) sin(z0) ],  z0 = 90.833 deg
    """
    return acosd((sind(decl) - sind(lat) * cosd(z0)) / (cosd(lat) * sind(z0)))


def sunset_azimuth(lat, decl, z0=90.833):
    return 360.0 - sunrise_azimuth(lat, decl, z0)


# ======================================================================
# verification against NOAA's own code, run under node
# ======================================================================
def _build_jobs(path, lat=SITE["lat"], lon=SITE["lon"], tz=SITE["tz"], year=2025):
    """Every day of `year`, sampled every 3 h of IST -- the benchmark grid."""
    jobs = []
    d0 = date(year, 1, 1)
    for doy in range(1, days_in_year(year) + 1):
        dt = d0 + timedelta(days=doy - 1)
        for hh in range(0, 24, 3):
            jobs.append({"id": f"{dt.year}-{dt.month:02d}-{dt.day:02d}T{hh:02d}",
                         "year": dt.year, "month": dt.month, "day": dt.day,
                         "lat": lat, "lon": lon, "tz": tz, "minutes_local": hh * 60.0})
    with open(path, "w") as f:
        json.dump(jobs, f)
    return jobs


def verify(workdir):
    jobs_path = os.path.join(workdir, "jobs.json")
    out_path = os.path.join(workdir, "oracle.json")
    jobs = _build_jobs(jobs_path)
    if not os.path.exists(out_path):
        subprocess.run(["node", os.path.join(workdir, "noaa_oracle.js"), jobs_path, out_path], check=True)
    oracle = {r["id"]: r for r in json.load(open(out_path))}

    import statistics

    def rms(v):  return math.sqrt(sum(x * x for x in v) / len(v)) if v else 0.0

    st = {k: [] for k in ("full_az", "pdf_az", "simple_az", "pdf_decl", "simple_decl", "pdf_eot", "simple_eot")}
    worst_doy = 0
    n_day = 0
    sampled = []
    for j in jobs:
        o = oracle[j["id"]]
        y, m, d = j["year"], j["month"], j["day"]
        doy = day_of_year(y, m, d)
        worst_doy = max(worst_doy, abs(doy - o["doy"]))
        # Azimuth is only meaningful/well conditioned for 5 deg < altitude < 85 deg:
        # below the horizon the Sun is heading for the nadir, and within ~1 deg of the
        # zenith every azimuth meets, so the azimuth loses all precision there.
        if o["el"] <= 5.0 or o["el"] >= 85.0:
            continue
        n_day += 1
        az1, _, _, decl1, eqt1 = noaa_pdf_azel(y, m, d, j["minutes_local"], j["lat"], j["lon"], j["tz"])
        az2, _, _, decl2, eqt2, _ = noaa_full_azel(julian_centuries(jd_from_ymd(y, m, d)),
                                                   j["minutes_local"], j["lat"], j["lon"], j["tz"])
        decl0, eot0 = simple_decl(doy), simple_eot(doy)
        H = hour_angle(clock_to_solar_time(j["minutes_local"] / 60.0, doy, j["lon"]))
        az0 = azimuth_general(j["lat"], decl0, H)
        st["full_az"].append(abs(wrap180(o["az"] - az2)))
        st["pdf_az"].append(abs(wrap180(az1 - o["az"])))
        st["simple_az"].append(abs(wrap180(az0 - o["az"])))
        st["pdf_decl"].append(abs(decl1 - o["decl"]))
        st["simple_decl"].append(abs(decl0 - o["decl"]))
        st["pdf_eot"].append(abs(eqt1 - o["eqtime"]))
        st["simple_eot"].append(abs(eot0 - o["eqtime"]))
        sampled.append({"doy": doy, "hour": j["minutes_local"] / 60.0, "el": o["el"], "az": o["az"],
                        "err_pdf": wrap180(az1 - o["az"]), "err_simple": wrap180(az0 - o["az"])})

    # ---- sunrise-azimuth benchmark over the whole year -------------------
    daily = []
    d0 = date(2025, 1, 1)
    for doy in range(1, 366):
        dt = d0 + timedelta(days=doy - 1)
        o = oracle[f"{dt.year}-{dt.month:02d}-{dt.day:02d}T06"]
        rise_min_local = o["rise_time_local"]
        a_noaa = o["rise_az"]
        a_simple = sunrise_azimuth(SITE["lat"], simple_decl(doy))
        rise_hr = (rise_min_local) / 60.0
        g = noaa_pdf_gamma(doy, rise_hr)
        a_pdf = sunrise_azimuth(SITE["lat"], noaa_pdf_decl(g))
        eot_pdf = noaa_pdf_eqtime(g)
        t_simple = solar_time_to_clock(12.0 - sunrise_hour_angle(SITE["lat"], simple_decl(doy)) / 15.0,
                                       doy, SITE["lon"]) * 60.0
        t_pdf = (12.0 - sunrise_hour_angle(SITE["lat"], noaa_pdf_decl(g)) / 15.0
                 - eot_pdf / 60.0 - (SITE["lon"] - SITE["lon_std"]) / 15.0) * 60.0
        daily.append({"date": dt.isoformat(), "doy": doy, "a_noaa": a_noaa, "a_simple": a_simple,
                      "a_pdf": a_pdf, "err_simple": a_simple - a_noaa, "err_pdf": a_pdf - a_noaa,
                      "rise_noaa": rise_min_local, "err_t_simple": t_simple - rise_min_local,
                      "err_t_pdf": t_pdf - rise_min_local})

    err_s = [abs(r["err_simple"]) for r in daily]
    err_p = [abs(r["err_pdf"]) for r in daily]
    rep = {
        "jobs": len(jobs), "daylight_samples": n_day, "worst_doy": worst_doy,
        "transcription_max_az": max(st["full_az"]),
        "pdf_az_max": max(st["pdf_az"]), "pdf_az_rms": rms(st["pdf_az"]),
        "simple_az_max": max(st["simple_az"]), "simple_az_rms": rms(st["simple_az"]),
        "pdf_decl_max": max(st["pdf_decl"]), "simple_decl_max": max(st["simple_decl"]),
        "pdf_eot_max": max(st["pdf_eot"]), "simple_eot_max": max(st["simple_eot"]),
        "rise_az_err_simple_max": max(err_s), "rise_az_err_simple_rms": rms(err_s),
        "rise_az_err_pdf_max": max(err_p), "rise_az_err_pdf_rms": rms(err_p),
        "rise_time_err_simple_max_min": max(abs(r["err_t_simple"]) for r in daily),
        "rise_time_err_pdf_max_min": max(abs(r["err_t_pdf"]) for r in daily),
        "annual_az_min": min(r["a_noaa"] for r in daily), "annual_az_max": max(r["a_noaa"] for r in daily),
    }
    with open(os.path.join(workdir, "verification.json"), "w") as f:
        json.dump({"summary": rep, "daily": daily, "sampled": sampled}, f, indent=1)

    print(f"benchmark grid : {len(jobs)} samples ({n_day} with 5 deg < altitude < 85 deg)")
    print(f"day-of-year    : Python vs NOAA, max |dN| = {worst_doy} day")
    print("--- Level 2  (full NOAA algorithm, transcription check) -----------")
    print(f"    max |dA| = {rep['transcription_max_az']:.6f} deg  (machine-identical to NOAA's own code)")
    print("--- Level 1  (NOAA published sheet, Sun up) -----------------------")
    print(f"    azimuth  max {rep['pdf_az_max']:.3f} deg, rms {rep['pdf_az_rms']:.3f} deg")
    print(f"    declination max {rep['pdf_decl_max']:.4f} deg | EoT max {rep['pdf_eot_max']:.3f} min")
    print("--- Level 0  (Class-12 model, Sun up) -----------------------------")
    print(f"    azimuth  max {rep['simple_az_max']:.3f} deg, rms {rep['simple_az_rms']:.3f} deg")
    print(f"    declination max {rep['simple_decl_max']:.4f} deg | EoT max {rep['simple_eot_max']:.3f} min")
    print("--- sunrise azimuth, all 365 days of 2025 (headline result) -------")
    print(f"    Class-12 formula : max {rep['rise_az_err_simple_max']:.3f} deg, rms {rep['rise_az_err_simple_rms']:.3f} deg")
    print(f"    NOAA sheet       : max {rep['rise_az_err_pdf_max']:.3f} deg, rms {rep['rise_az_err_pdf_rms']:.3f} deg")
    print(f"    sunrise time     : Class-12 max {rep['rise_time_err_simple_max_min']:.2f} min,"
          f" NOAA sheet max {rep['rise_time_err_pdf_max_min']:.2f} min")
    print(f"    NOAA sunrise azimuth range over the year: {rep['annual_az_min']:.2f}..{rep['annual_az_max']:.2f} deg")
    print("    -> written to work/verification.json")


# ======================================================================
# tables
# ======================================================================
def _fmt_time(minutes):
    minutes = minutes % 1440.0
    h = int(minutes // 60)
    m = int(round(minutes - h * 60))
    if m == 60:
        h, m = (h + 1) % 24, 0
    return f"{h:02d}:{m:02d}"


def annual_table(y, out_csv, lat=SITE["lat"], lon=SITE["lon"], tz=SITE["tz"]):
    rows = []
    for m in range(1, 13):
        for d in (1, 8, 15, 22):
            rt, ra, decl = noaa_full_rise_set(True, y, m, d, lat, lon, tz)
            st, sa, _ = noaa_full_rise_set(False, y, m, d, lat, lon, tz)
            n = day_of_year(y, m, d)
            rows.append({
                "date": f"{y}-{m:02d}-{d:02d}", "N": n,
                "decl_deg": round(decl, 4),
                "decl_simple_deg": round(simple_decl(n), 4),
                "eot_min": round(noaa_full_eqtime(julian_centuries(jd_from_ymd(y, m, d))), 3),
                "eot_simple_min": round(simple_eot(n), 3),
                "sunrise_IST": _fmt_time(rt), "sunrise_az_deg": round(ra, 3),
                "sunset_IST": _fmt_time(st), "sunset_az_deg": round(sa, 3),
                "solar_noon_IST": _fmt_time(noaa_full_solar_noon(y, m, d, lon, tz)),
            })
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


def region_table(out_csv, lat=SITE["lat"], lon=SITE["lon"], tz=SITE["tz"], radius_km=REGION_RADIUS_KM):
    """How much does sunrise azimuth/time change across the 22 km disc?"""
    dlat, dlon = region_bounds(lat, lon, radius_km)
    corners = {
        "centre": (lat, lon),
        "N edge": (lat + dlat, lon),
        "S edge": (lat - dlat, lon),
        "E edge": (lat, lon + dlon),
        "W edge": (lat, lon - dlon),
        "NE corner": (lat + dlat, lon + dlon),
        "SW corner": (lat - dlat, lon - dlon),
    }
    rows = []
    for (y, m, d, label) in [(2025, 6, 21, "Jun solstice"), (2025, 12, 21, "Dec solstice"),
                             (2025, 3, 20, "Mar equinox"), (2025, 9, 22, "Sep equinox")]:
        for name, (la, lo) in corners.items():
            rt, ra, _ = noaa_full_rise_set(True, y, m, d, la, lo, tz)
            rows.append({"event": label, "point": name, "lat": round(la, 6), "lon": round(lo, 6),
                         "sunrise_IST": _fmt_time(rt), "sunrise_min": round(rt, 3),
                         "sunrise_az_deg": round(ra, 4)})
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return rows


def bound_table(out_csv):
    dlat, dlon = region_bounds(SITE["lat"], SITE["lon"])
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["bound", "lat", "lon", "km from centre (N-S)", "km from centre (E-W)"])
        w.writerow(["north", f"{SITE['lat']+dlat:.6f}", f"{SITE['lon']:.6f}", 22.0, 0.0])
        w.writerow(["south", f"{SITE['lat']-dlat:.6f}", f"{SITE['lon']:.6f}", -22.0, 0.0])
        w.writerow(["east", f"{SITE['lat']:.6f}", f"{SITE['lon']+dlon:.6f}", 0.0, 22.0])
        w.writerow(["west", f"{SITE['lat']:.6f}", f"{SITE['lon']-dlon:.6f}", 0.0, -22.0])
    return dlat, dlon


# ======================================================================
# demo
# ======================================================================
def demo():
    lat, lon = SITE["lat"], SITE["lon"]
    y, m, d = 2025, 12, 21
    n = day_of_year(y, m, d)
    print("=" * 74)
    print(f"WORKED EXAMPLE -- winter solstice {y}-{m:02d}-{d:02d} at {lat:.4f} N {lon:.4f} E")
    print("=" * 74)
    decl0, eot0 = simple_decl(n), simple_eot(n)
    print(f"day number N                        = {n}")
    print(f"declination (simple)                = {decl0:.3f} deg")
    print(f"equation of time (simple)           = {eot0:.3f} min")
    z0 = 90.833
    H0 = sunrise_hour_angle(lat, decl0, z0)
    A = sunrise_azimuth(lat, decl0, z0)
    print(f"sunrise hour angle H0               = {H0:.4f} deg = {H0/15:.4f} h")
    print(f"sunrise azimuth (simple + z0)       = {A:.3f} deg   -> compass {compass(A)}")
    rt, ra, decl2 = noaa_full_rise_set(True, y, m, d, lat, lon, SITE["tz"])
    print(f"sunrise azimuth (NOAA)              = {ra:.3f} deg   -> compass {compass(ra)}")
    print(f"sunrise clock time (simple)         = {_fmt_time(solar_time_to_clock(12 - H0/15, n, lon)*60)} IST")
    print(f"sunrise clock time (NOAA)           = {_fmt_time(rt)} IST")
    print(f"declination (NOAA)                  = {decl2:.3f} deg")
    print(f"solar noon (NOAA)                   = {_fmt_time(noaa_full_solar_noon(y, m, d, lon, SITE['tz']))} IST")
    print("-" * 74)
    print("hourly azimuth on the same day (general equation):")
    for hh in (6.5, 8, 10, 12, 14, 16, 17.5):
        H = hour_angle(clock_to_solar_time(hh, n, lon))
        print(f"   {hh:5.1f} h IST -> A = {azimuth_general(lat, decl0, H):7.3f} deg"
              f"   (NOAA {noaa_full_azel(julian_centuries(jd_from_ymd(y,m,d)), hh*60, lat, lon, SITE['tz'])[0]:7.3f} deg)")


def compass(a):
    names = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return names[int((a % 360) / 22.5 + 0.5) % 16]


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "demo"
    work = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work")
    os.makedirs(work, exist_ok=True)
    if cmd == "verify":
        verify(work)
    elif cmd == "tables":
        os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"), exist_ok=True)
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        annual_table(2025, os.path.join(d, "annual_2025.csv"))
        region_table(os.path.join(d, "region_22km.csv"))
        bound_table(os.path.join(d, "region_bounds.csv"))
        print("tables written to", d)
    else:
        demo()
