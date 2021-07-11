import os, sys
import numpy as np

from astropy.table import Table

rad2deg = 180 / np.pi

def random_positions(ra_range, dec_range, Ncoord):
    ra1, ra2 = np.sort(ra_range)
    dec1, dec2 = np.sort(dec_range)

    RA = np.random.uniform(low=ra1, high=ra2, size=Ncoord)
    x = np.random.uniform(low=np.sin(dec1/rad2deg),\
                          high=np.sin(dec2/rad2deg),\
                          size=Ncoord)
    Dec  = np.arcsin(x) * rad2deg

    return RA, Dec

def cartesian_to_polar(x, y):
    z = x + y*1j
    r = np.abs(z)
    theta = np.angle(z) * rad2deg

    return r, theta


