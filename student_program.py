# sun_azimuth.py -- "Where is the Sun right now?"  (Class XII, Python 3)
# Site: 23.1698 N, 82.3603 E (IST = UTC+5:30).  Angles in degrees.
from math import sin, cos, tan, asin, acos, atan2, radians, degrees

LAT, LON, TZ, LON_STD = 23.1698, 82.3603, 5.5, 82.5      # site and IST meridian

def day_number(d, m, y):                       # 1 January = day 1
    cum = [0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    leap = (y % 4 == 0 and y % 100 != 0) or y % 400 == 0
    return cum[m] + d + (1 if (m > 2 and leap) else 0)

def declination(n):                            # Eq. (8):  Cooper (1969)
    return 23.44 * sin(radians(360.0 / 365.0 * (n - 81)))

def equation_of_time(n):                       # Eq. (9):  minutes
    B = radians(360.0 * (n - 81) / 364.0)
    return 9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B)

def sun_azimuth(clock_hours, n):
    """Azimuth of the Sun (degrees clockwise from true North) and its altitude."""
    solar_time = clock_hours + equation_of_time(n) / 60.0 + (LON - LON_STD) / 15.0
    H = 15.0 * (solar_time - 12.0)             # Eq. (10): hour angle, degrees
    d = radians(declination(n))                # declination
    p = radians(LAT)                           # latitude
    h = radians(H)
    A = 180.0 + degrees(atan2(sin(h), cos(h) * sin(p) - tan(d) * cos(p)))   # Eq. (13)
    z = degrees(acos(sin(p) * sin(d) + cos(p) * cos(d) * cos(h)))           # zenith angle
    return A % 360.0, 90.0 - z

def sunrise(clock_date=None, n=None):
    """Sunrise azimuth, Eq. (12), and sunrise clock time, Eq. (14)+EoT."""
    n = n or day_number(clock_date[0], clock_date[1], clock_date[2])
    d, p, z0 = radians(declination(n)), radians(LAT), radians(90.833)
    A = degrees(acos((sin(d) - sin(p) * cos(z0)) / (cos(p) * sin(z0))))
    H0 = degrees(acos(cos(z0) / (cos(p) * cos(d)) - tan(p) * tan(d)))
    t_solar = 12.0 - H0 / 15.0                                     # sunrise, solar time
    t_clock = t_solar - equation_of_time(n) / 60.0 - (LON - LON_STD) / 15.0
    return A, t_clock

def hhmm(x):
    return "%02d:%02d" % (int(x), round(60 * (x - int(x))))

n = day_number(21, 12, 2025)
print("21 Dec 2025, site 23.1698 N 82.3603 E")
print("  declination      = %+.3f deg" % declination(n))
print("  equation of time = %+.2f min" % equation_of_time(n))
for t in (6, 8, 10, 12, 14, 16, 18):
    A, alt = sun_azimuth(t, n)
    print("  %02d:00 IST  azimuth %7.3f deg   altitude %6.2f deg" % (t, A, alt))
A, tr = sunrise((21, 12, 2025))
print("  sunrise: %s IST, azimuth %.3f deg (i.e. %.3f deg south of due east)" % (hhmm(tr), A, A - 90))
