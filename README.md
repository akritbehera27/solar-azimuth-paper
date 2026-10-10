# Where Does the Sun Rise? — code and data

Supporting material for the paper

> **Where Does the Sun Rise?**
> A general equation for the Sun's azimuth — built on NOAA's verified solar-position equations,
> in the views of a Class-XII student, and applied to a 22 km region around
> 23.169753993380713° N, 82.36029982215409° E.

**Author:** Akrit Behera (Kendriya Vidyalaya, Class XII)
**ORCID:** https://orcid.org/0009-0004-6530-1901
**Contact:** akritbehera27@gmail.com

---

## Files

| File | What it is | Size |
|---|---|---|
| `solar_azimuth.py` | The complete computing engine. Three levels of approximation (NOAA full Fourier series → NOAA's printed-sheet constants → a Class-XII-level model), plus sunrise/sunset azimuth routines, a worked example, table generators, and a self-check that reproduces every residual in the paper. Run with `python3 solar_azimuth.py verify`, `python3 solar_azimuth.py tables`, `python3 solar_azimuth.py demo`. | ~25 KB |
| `student_program.py` | The ~60-line, dependency-free Python program printed in Appendix A of the paper. Pure textbook-level code; suitable for a CBSE Class-XII computer-science project. | ~2.6 KB |
| `annual_2025.csv` | For each of 47 representative dates of 2025: calendar date, solar declination, equation of time, sunrise azimuth (degrees), sunrise IST, sunset azimuth, sunset IST, solar noon IST, day length. The table printed as Table 1 in the paper. | ~3.9 KB |
| `region_22km.csv` | The same quantities evaluated at the centre of the 22 km disc and at its four cardinal edges (N/S/E/W), demonstrating the invariance result. Table 4 in the paper. | ~1.8 KB |
| `region_bounds.csv` | The latitude and longitude of the north, south, east and west edges of the 22 km disc around the study site, in degrees. | 201 B |
| `verification.json` | The raw residual data behind every accuracy claim in the paper: 2920-sample grid of NOAA full vs. simplified vs. textbook-level azimuths, declinations, EoTs, and sunrise times. | ~332 KB |

---

## Reproducing the paper's numbers

```bash
# 1. Cross-check every model against NOAA's own algorithm
python3 solar_azimuth.py verify

# 2. Rebuild the CSVs used for the tables
python3 solar_azimuth.py tables

# 3. Run the worked example (interior of the Sun, 21 December 2025)
python3 solar_azimuth.py demo
```

All three commands are deterministic and require only a standard Python 3 installation
(no third-party packages).

---

## Source of the equations

The reference equations are those of NOAA's *General Solar Position Calculations*
(NOAA GML, n.d.-a) and of the NOAA Solar Calculator (NOAA GML, n.d.-b), available at

- https://gml.noaa.gov/grad/solcalc/
- https://gml.noaa.gov/grad/solcalc/solareqns.PDF

NOAA's own published JavaScript implementation (`work/noaa_main.js` in the paper's
repository) was downloaded unmodified and executed for 2920 instants as the numerical
reference. The simplified model uses only the leading harmonics of NOAA's Fourier series,
as is standard in the solar-engineering literature (Cooper, 1969).

---

## License

This code and data are released under the Creative Commons Attribution 4.0 International
licence (CC BY 4.0). You may share and adapt the material for any purpose, including
commercially, as long as you credit the author.

https://creativecommons.org/licenses/by/4.0/

