"""Sensitivity information from Khostovan+2020"""
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from icecream import ic
from matplotlib import pyplot as plt
import numpy as np


z = 0.47
Da = cosmo.angular_diameter_distance(z)
# required luminosity corrected by [NII] contamination
LHa = 1.2e40 * u.erg / u.s / 10**0.11
ic(LHa)
FHa = LHa / (4*np.pi*Da**2)
FHa = FHa.to(u.erg/u.s/u.cm**2)
ic(FHa)

texp_cosmos = 47.25*u.h
limit_cosmos = 8.2e-18

texp = np.linspace(5*u.min, texp_cosmos, 100)
limit = limit_cosmos * (texp_cosmos/texp)**0.5

j = np.argmin(np.abs(limit - FHa.value))
ic(texp[j], limit[j])

fig, ax = plt.subplots()
ax.plot(texp, limit, 'C0-')
ax.plot(texp_cosmos, limit_cosmos, 'C1o')
ax.plot(texp[j], limit[j], 'C2*')
ax.set(yscale='log', xlabel='Exposure time (h)',
       ylabel=r'$5\sigma$ limit (erg cm$^{-2}$ s$^{-1}$)')
plt.show()