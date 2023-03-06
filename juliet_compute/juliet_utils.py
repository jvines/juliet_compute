"""Utilities specific to Juliet."""
import pickle

from .globals import *
from .planetary import *


def read_juliet_posteriors(file):
    """Read a posterior pickle file from Juliet.

    Parameters
    ----------
    file: str
        The path to the Juliet posterior file.

    Returns
    -------
    res: dict
        A dictionary with the relevant parameter samples.
    """
    out = dict()
    with open(file, 'rb') as jar:
        res = pickle.load(jar)
    try:
        out['pl'], out['pu'] = res['pl'], res['pu']
    except KeyError:
        pass
    posts = res['posterior_samples']
    for k in posts.keys():
        for param in params_of_interest:
            if f'{param}_p' in k:
                out[k] = posts[k]
                break
        if 'q1' in k or 'q2' in k or 'rho' == k:
            out[k] = posts[k]
    return out


def get_number_of_planets(post):
    """Get total number of planets from the posteriors."""
    # TODO FIX
    planets = list()
    for param in post.keys():
        planet = param.split('_')[-1]
        if planet not in planets:
            planets.append(planet)
    return len(planets)


def get_bp(post, nplanets, pl, pu):
    """Get impact parameter and radius ratio from r1 and r2."""
    for npl in range(1, nplanets + 1):
        if f'r1_p{npl}' in post:
            b, p = obtain_bp(post[f'r1_p{npl}'], post[f'r2_p{npl}'], pl, pu)
            post[f'b_p{npl}'] = b
            post[f'p_p{npl}'] = p
    return True


def get_radius(post, nplanets, unit):
    """Get radius of planets from the posteriors."""
    # If no stellar radius provided, abort.
    if post['Rs'][0] == -1:
        return False
    for npl in range(1, nplanets + 1):
        if f'p_p{npl}' in post:
            r = radius(post['Rs'], post[f'p_p{npl}'], unit)
            post[f'rp_p{npl}'] = r
    return True


def get_semimajor_axis(post, nplanets):
    """Calculate the scaled and non scaled semimajor axis in AU."""
    # If no stellar radius provided, abort.
    if post['Rs'][0] == -1:
        return False
    if 'rho' in post.keys():  # if rho parametrization was used.
        rho = post['rho']
        for npl in range(1, nplanets + 1):
            # Calculate scaled semimajor axis.
            a = rho_to_a(rho, post[f'P_p{npl}'])
            post[f'a_p{npl}'] = a  # Save the scaled result.
            # Transform to AU.
            aau = a * post['Rs'] * solR2AU
            post[f'aau_p{npl}'] = aau  # Save the AU result.
    else:  # rho wasn't used so we already have a_pN
        for npl in range(1, nplanets + 1):
            post[f'aau_p{npl}'] = post[f'a_p{npl}'] * post['Rs'] * solR2AU
    return True


def get_inclination(post, nplanets):
    """Get the system inclination angles in radians."""
    for npl in range(1, nplanets + 1):
        if f'b_p{npl}' in post:  # if impact parameter for a planet.
            i = inclination(post[f'b_p{npl}'], post[f'a_p{npl}'])
            post[f'inc_p{npl}'] = i
    return True


def get_mass(post, unit, nplanets):
    """Calculate mass."""
    # If no stellar mass provided, abort.
    if post['Ms'][0] == -1:
        return False
    planets_with_amplitudes = list()
    full_mass = {f'p{n}': False for n in range(1, nplanets + 1)}
    # Get which planets have amplitude and impact param (and thus inclination).
    for param in post:
        pln = param.split('_')[-1]
        if 'K_p' in param:
            planets_with_amplitudes.append(pln)
        if 'inc_p' in param:
            full_mass[pln] = True

    if not len(planets_with_amplitudes):  # If no amplitudes are found.
        return False
    # For each planet calculate mass and append to post.
    for pln in full_mass.keys():
        if full_mass[pln]:
            marr = mass(
                post['Ms'],
                post[f'K_p{pln}'],
                post[f'P_p{pln}'],
                post[f'ecc_p{pln}'],
                post[f'inc_p{pln}'],
                unit=unit
            )
            post[f'mp_{pln}'] = marr
        else:
            marr = msini(
                post['Ms'],
                post[f'K_p{pln}'],
                post[f'P_p{pln}'],
                post[f'ecc_p{pln}'],
                unit=unit
            )
            post[f'msinip_{pln}'] = marr
    return True


def get_teq(post, nplanets, ab):
    """Calculate equilibrium temperature."""
    # TODO
    # If no stellar radius provided, abort.
    if post['Rs'][0] == -1:
        return False
    # If no effective temperature provided, abort.
    if post['Teff'][0] == -1:
        return False
    for npl in range(1, nplanets + 1):
        pass
    return True


def get_transit_duration(post, nplanets):
    """Calculate transit duration."""
    # TODO
    # If no stellar radius provided, abort.
    if post['Rs'][0] == -1:
        return False
    for npl in range(1, nplanets + 1):
        pass
    return True
