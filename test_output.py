import os, sys

import numpy as np
from astropy.table import Table

from lensing import lens_one_object
import utils

def main():
    tbl = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/highz/allr.fits')

    for index in range(len(tbl)):
        # quasar coordinate
        qra = tbl['qra'][index]
        qdec = tbl['qdec'][index]

        # galaxy coordinate
        gra = tbl['ra'][index]
        gdec = tbl['dec'][index]

        # galaxy AR and PA
        ar = tbl['axisratio'][index]
        pa = tbl['pa_lens'][index]

        # galaxy redshift
        zl = tbl['redshift_true'][index]

        # quasar redshift
        zs = tbl['z'][index]

        # galaxy velosity dispersion
        sigma = tbl['sigma'][index]

        # externel convergence and shear
        convergence = tbl['convergence'][index]
        shear_1 = tbl['shear_1'][index]
        shear_2 = tbl['shear_2'][index]

        # convert cartesian to polar
        shear, pa_shear = utils.cartesian_to_polar(shear_1, shear_2)

        # only include targets with lens redshift lower than source
#        if zs > zl+0.01:
        if True:
            image_x, image_y, image_mu, image_dt = \
                lens_one_object(qra, qdec, gra, gdec,\
                                zs, zl, sigma, ar, pa,\
                               convergence, shear, pa_shear)

        else:
            image_x, image_y, image_mu, image_dt = \
                    [np.array([-999.])] * 4


        print(image_x)
        print(tbl['image_x'][index])

        print(image_y)
        print(tbl['image_y'][index])

        print(image_mu)
        print(tbl['image_mu'][index])

        flag = input('Continue?')
        if flag in 'Yy':
            continue

        else:
            break

if __name__=='__main__':
    main()
