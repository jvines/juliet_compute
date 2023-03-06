"""Global parameters."""
import astropy.constants as const
from astropy import units
from numpy import log, pi

# USEFUL CONSTANTS
r2d = 180 / pi
d2r = 1 / r2d
solR2AU = units.solRad.to(units.AU)
solR2J = units.solRad.to(units.jupiterRad)
solR2E = units.solRad.to(units.earthRad)
G = const.G.to(units.cm ** 3 / units.g / units.d ** 2)
# Transform from Mj/Rj^3 to g/cm^3
jupden2gcc = (1 * units.jupiterMass / units.jupiterRad ** 3).cgs.value

# Juliet specific
params_of_interest = [
    'P', 'K', 'p', 'b', 't0', 'r1', 'r2', 'ecc', 'omega',
    'a', 'esinomega', 'ecosomega', 'sesinomega', 'secosomega'
]
