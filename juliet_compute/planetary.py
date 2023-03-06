"""This module holds planetary and orbital functions."""
import astropy.units as u
import numpy as np

from .error import UnitError
from .globals import G, d2r, solR2J, solR2E


# MASS FUNCTIONS #


def msini(Ms, K, P, e, unit='jupiter'):
    """Compute M sini from stellar mass, amplitude, period and eccentricity.

    Parameters
    ----------
    Ms: array_like
        An array of floats containing the stellar mass in solar radii.
    K: array_like
        An array of floats containing the amplitude in m/s.
    P: array_like
        An array of floats containing the period in days.
    e: array_like
        An array of floats containing the eccentricity.
    unit: str, optional
        A string indicating the output units. Can be `jupiter` or `earth`.

    Returns
    -------
    ret: array_like
        An array of floats containing the minimum mass.
    """
    pp = (P * u.day) ** (1 / 3)
    ee = np.sqrt(1 - e ** 2)
    mm = (Ms * u.solMass) ** (2 / 3)
    msini = K * (u.m / u.s) * ee * pp * mm / (2 * np.pi * G) ** (1 / 3)
    if unit == 'jupiter':
        ret = msini.to(u.jupiterMass).value
    elif unit == 'earth':
        ret = msini.to(u.earthMass).value
    else:
        UnitError('Mass').__raise__()
    return ret


def mass(Ms, K, P, e, i, unit='jupiter'):
    """Compute the mass from stellar mass, amplitude, period and eccentricity.

    Parameters
    ----------
    Ms: array_like
        An array of floats containing the stellar mass in solar radii.
    K: array_like
        An array of floats containing the amplitude in m/s.
    P: array_like
        An array of floats containing the period in days.
    e: array_like
        An array of floats containing the eccentricity.
    i: array_like
        An array of floats containing the system inclination in degrees.
    unit: str, optional
        A string indicating the output units. Can be `jupiter` or `earth`.

    Returns
    -------
    ret: array_like
        An array of floats containing the mass.
    """
    pp = (P * u.day) ** (1 / 3)
    ee = np.sqrt(1 - e ** 2)
    mm = (Ms * u.solMass) ** (2 / 3)
    sini = np.sin(d2r * i)
    m = K * (u.m / u.s) * ee * pp * mm / (2 * np.pi * G) ** (1 / 3) / sini
    if unit == 'jupiter':
        ret = m.to(u.jupiterMass).value
    elif unit == 'earth':
        ret = m.to(u.earthMass).value
    else:
        UnitError('Mass').__raise__()
    return ret


# RADIUS FUNCTIONS #

def radius(Rs, p, unit='jupiter'):
    """Compute planetary radius from stellar radius and radius ratio.

    Parameters
    ----------
    Rs: array_like
        An array of floats corresponding to the stellar radius.
    p: array_like
        An array of floats corresponding to the radius ratio.
    unit: str, optional
        What unit to use for the planetary radius. Can be `jupiter` or `earth`

    Returns
    -------
    rp: array_like
        An array of floats corresponding to the planetary radius.
    """
    if unit == 'jupiter':
        rp = p * Rs * solR2J
    elif unit == 'earth':
        rp = p * Rs * solR2E
    else:
        UnitError('Radius').__raise__()
    return rp


def obtain_bp(r1, r2, pl, pu):
    """Transform from pu, pl and Ar to impact parameter and radius ratio.

    Parameters
    ----------
    r1: array_like
        An array of floats containing the r1 parametrization.
    r2: array_like
        An array of floats containing the r2 parametrization.
    pl: array_like
        An array of floats containing the lower radius ratio limit.
    pu: array_like
        An array of floats containing the upper radius ratio limit.

    Returns
    -------
    b: array_like
        The impact parameter.
    p: array_like
        The radius ratio.
    """
    ar = (pu - pl) / (2. + pl + pu)
    if r1 > ar:
        b = (1 + pl) * (1 + (r1 - 1) / (1 - ar))
        p = (1 - r2) * pl + r2 * pu
    else:
        b = (1 + pl) + np.sqrt(r1 / ar) * r2 * (pu - pl)
        p = pu + (pl - pu) * np.sqrt(r1 / ar) * (1 - r2)
    return b, p


def density(mass, radius, unit='jupiter'):
    aunit = u.jupiterMass / u.jupiterRad ** 3
    if unit == 'earth':
        aunit = u.earthMass / u.earthRad ** 3
    den = mass / (4 / 3 * np.pi * radius ** 3) * aunit.to(u.g / u.cm ** 3)
    return den

# EFFECTIVE TEMPERATURE FUNCTIONS #


def teq(Teff, Rs, Ab, a):
    """Compute the equilibrium temperature.

    Parameters
    ----------
    Teff: array_like
        An array of floats containing the stellar effective temperature
        in Kelvin.
    Rs: array_like
        An array of floats containing the stellar radius in solar Radii.
    Ab: float
        The albedo.
    a: array_like
        An array of floats containing the orbital semimajor axis in solar Radii.

    Returns
    -------
    teq: array_like
        The equilibrium temperature in Kelvin.
    """
    return Teff * np.sqrt(Rs / (2 * a)) * (1 - Ab) ** (1 / 4)


# ORBITAL TRANSFORMATIONS #

def get_ew(s, c, parametrization='sesinw'):
    if parametrization == 'sesinw':
        e = s ** 2 + c ** 2
    elif parametrization == 'esinw':
        e = (s ** 2 + c ** 2) ** .5
    w = np.arccos(c / e ** .5)
    return e, w


def rho_to_a(rho, P):
    """Compute scaled semimajor axis.

    Parameters
    ----------
    rho: array_like
        An array of floats containing the stellar density in g/cm3.
    P: array_like
        An array of floats containing the orbital period in days.

    Returns
    -------
    a: array_like
        An array containing the semimajor axis (a/Rs).
    """
    return (rho * G.value * P ** 2 / (3 * np.pi)) ** (1 / 3)


def true_anomaly(w):
    """Compute the true anomaly.

    Parameters
    ----------
    w: array_like
        An array of floats containing the longitude of periastron or
        argument of periastron passage in radians.

    Returns
    -------
    f: array_like
        An array of floats containing the true anomaly in radians.
    """
    return np.pi / 2 - w


def eccentric_anomaly(e, f):
    """Compute the eccentric anomaly.

    Parameters
    ----------
    e: array_like
        An array of floats containing the eccentricity.
    f: array_like
        An array of floats containing the true anomaly in radians.

    Returns
    -------
    E: array_like
        An array of floats containing the eccentric anomaly in radians.
    """
    return 2 * np.arctan(np.sqrt(1 - e) / (1 + 2) * np.tan(f / 2))


def mean_anomaly(e, E):
    """Compute the mean anomaly.

    Parameters
    ----------
    e: array_like
        An array of floats containing the eccentricity.
    E: array_like
        An array of floats containing the eccentric anomaly in radians.

    Returns
    -------
    M: array_like
        An array of floats containing the mean anomaly in radians.
    """
    return E - e * np.sin(E)


def transit_duration(P, Rs, Rp, a, e, w):
    """Compute transit duration in hours.

    Parameters
    ----------
    P: array_like
        An array of floats containing the orbital period in days.
    Rs: array_like
        An array of floats containing the stellar radius in solar Radii.
    Rp: array_like
        An array of floats containing the planetary radius in solar Radii.
    a: array_like
        An array of floats containing the semimajor axis in solar Radii.
    e: array_like
        An array of floats containing the eccentricity.
    w: array_like
        An array of floats containing the longitude of periastron or
        argument of periastron passage in radians.

    Returns
    -------
    dur: array_like
        An array of floats containing the transit duration in hours.
    """
    f = true_anomaly(w)
    E = eccentric_anomaly(e, f)
    df = np.arcsin((Rs + Rp) / a / (1 - e * np.cos(E)))
    f_adj = f + df
    E_adj = eccentric_anomaly(e, f_adj)
    M = mean_anomaly(e, E_adj)
    dur = P * M / np.pi / 24  # divide by 24 transforms from days to hours.
    return dur


def inclination(b, a):
    """Calculate orbital inclination.

    Parameters
    ----------
    b: array_like
        An array of floats containing the impact parameter.
    a: array_like
        An array of floats containing the scaled semimajor axis.

    Returns
    -------
    i: array_like
        An array of floats containing the inclination in radians.
    """
    i = np.arccos(b / a)
    return i


# LD COEFFICIENT TRANSFORMATIONS #
# References: Kipping 2013, Espinoza E. & Jordan A 2016


def q_to_u(law, q1, q2):
    """Transform kipping coefficients into u coefficients."""
    if law == 'linear':
        u1 = q1
        u2 = q2
    elif law == 'quadratic':
        u1 = 2 * np.sqrt(q1) * q2
        u2 = np.sqrt(q1) * (1 - 2 * q2)
    elif law == 'square-root':
        u1 = np.sqrt(q1) * (1 - 2 * q2)
        u2 = 2 * np.sqrt(q1) * q2
    elif law == 'logarithmic':
        u1 = 1 - np.sqrt(q1) * q2
        u2 = 1 - np.sqrt(q1)
    return u1, u2


def u_to_q(law, u1, u2):
    """Transform u coefficients into kipping coefficients."""
    if law == 'quadratic':
        q1 = (u1 + u2) ** 2
        q2 = u1 / (2 * (u1 + u2))
    elif law == 'square-root':
        q1 = (u1 + u2) ** 2
        q2 = u2 / (2 * (u1 + u2))
    elif law == 'logarithmic':
        q1 = (1 - u2) ** 2
        q2 = (1 - u1) / (1 - u2)
    return q1, q2
