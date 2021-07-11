import os, sys
import numpy as np
from astropy.table import Table, vstack, hstack
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from collections import Counter

import glafic as sim
import utils as utils


defaultcosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def theta_E(sigma, zl, zs, cosmo=defaultcosmo):
    '''
    Parameters:
        sigma: the velocity dispersion of the deflector, in km s-1
        zl: the redshift of the deflector (lens)
        zs: the redshift of the source
        cosmo: the cosmology model (astropy.cosmology.FlatLambdaCDM class)

    Returns:
        the Einstein radius, in arcsecond
    '''
    Ds = cosmo.angular_diameter_distance(zs)
    Dl = cosmo.angular_diameter_distance(zl)
    Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)

    theta = Dls / Ds * 4 * np.pi * (sigma/3e5)** 2 * 206265
    return theta

def maxsep(xlist, ylist):
    '''
    Obtain the maximum separation between arbitrary positions in a list

    Parameters:
        xlist: array-like
            The x-coordinates of the positions

        ylist: array-like
            The y-coordinates of the positions

    returns: maxsep
        The maximum separation.
    '''
    if not len(xlist)==len(ylist):
        raise ValueError('The input x-coord and y-coord must have the same length.')

    if len(xlist)<=1:
        return -1

    sep = []
    for i in range(len(xlist)):
        for j in range(i, len(xlist)):
            x1 = xlist[i]
            y1 = ylist[i]
            x2 = xlist[j]
            y2 = ylist[j]
            sep.append(np.sqrt((x1-x2)**2+(y1-y2)**2))

    return np.max(sep)

def fill_with_default(arr, length, fill=-999.):
    '''
    Fill an array with default values if the array is shorter than the given length.
    '''

    addlength = length - len(arr)
    if addlength<1:
        return arr
    else:
        addarr = np.full(addlength, fill)
        return np.append(arr, addarr)

def match_pairs(qsotbl, galtbl, nn=1):
    '''
    Find close galaxy-quasar pairs.

    Parameters:
        qsotbl: astropy.table.Table class
            The table that contains the coordinates of the quasars.
            Coordinate keys are qra and qdec.

        galtbl: astropy.table.Table class
            The table that contains the coordinates of teh galaxies.
            Coordinate keys are ra and dec.

        nn: int
            This function will select the nn-th match.

    Return: selected_qsoidx, selected_galidx
        The indice of matched quasars and galaxies
    '''

    # Generate indice for the two tables
    qsoidx_all = np.arange(len(qsotbl))
    galidx_all = np.arange(len(galtbl))

    # the velocity dispersions
    sigma = galtbl['sigma']

    # this is the largest possible Einstein radius, in arcsec
    theta_max = 4 * np.pi * (sigma/3e5)** 2 * 206265

    # define the coordinate variables
    qsocoord = SkyCoord(ra=qsotbl['qra'], dec=qsotbl['qdec'], unit='deg')
    galcoord = SkyCoord(ra=galtbl['ra'], dec=galtbl['dec'], unit='deg')

    # matching the two catalogs.
    # idx: the indices into qsocoord that are the nn-th closest objects to each of the coordinates in galcoord
    # d2d: the on-sky distance 
    # d3d: dummy
    idx, d2d, d3d = galcoord.match_to_catalog_sky(qsocoord, nthneighbor=nn)

    # get matched pairs
    # we do not consider pairs that are 5 theta_E_max away
    matched = (d2d.to('arcsec').value < theta_max*5)

    # pick up the matched galaxies and quasars
    selected_galidx = galidx_all[matched]
    selected_qsoidx = idx[matched]

    return (selected_qsoidx, selected_galidx)

def save_crude_catalog(qsotbl, galtbl, seed=None):
    """
    Match the random quasar positions to the galaxy table.
    Specifically, for each galaxy, we find the closest quasar,
    and see if it might be lensed using an Einstein radius criteria.

    Parameters:
        qsotbl: astropy.table.Table
            The quasar table

        galtbl: astropy.table.Table
            The galaxy table

    Return: combined
        An astropy.table.Table variable, with each row a matched qso-gal pair that might generate a lens
    """

    if seed is not None:
        np.random.seed(seed)

    selected_qsoidx = []
    selected_galidx = []

    # start from the closest match, and loop until all matched distances are >5*theta_E

    nn = 1
    while 1:
        qi, gi = match_pairs(qsotbl, galtbl ,nn)
        if len(qi) < 1:
            break
        selected_qsoidx = selected_qsoidx + list(qi)
        selected_galidx = selected_galidx + list(gi)

        nn = nn + 1

    selected_qsotbl = qsotbl[:][selected_qsoidx]
    selected_galtbl = galtbl[:][selected_galidx]

    # combine info
    combined = hstack([selected_qsotbl,  selected_galtbl])
    combined['axisratio'] =\
            combined['size_minor_true'] / combined['size_true']

    # need to add random shears, convergences, PA
    shear_sig = np.log10(1+combined['z']) * 0.04024618 \
            - 0.00095052
    convergence_sig = np.log10(1+combined['z']) * 0.0566969 \
            - 0.00139201
    combined['pa_lens'] = np.random.uniform(0, 360, len(combined))
    combined['shear_1'] = np.random.normal(0, shear_sig, len(combined))
    combined['shear_2'] = np.random.normal(0, shear_sig, len(combined))
    combined['convergence'] =\
            np.random.normal(0, convergence_sig, len(combined))

    return combined

def lens_one_object(qra, qdec, gra, gdec, zs, zl, sigma, ar, pa,\
                   convergence, shear, pa_shear):

    '''
    Run GLAFIC to lens one system
    To test, compare the GLAFIC input file to the desired input
    '''

    # define coordinates of the background quasar, following the conversion in GLAFIC
    # I set the position of the foreground galaxy to be at (0,0)
    # x: ra-direction offset in arcsec
    # y: dec-direction offset in arcsec
    xqso = (gra - qra) * np.cos(gdec*np.pi / 180) * 3600
    yqso = (qdec - gdec) * 3600

    # define the lensing system, using an SIE lens
    lens_sie = {'Type':'sie', 'sigma':sigma, 'r_core':0.0,\
               'ellipticity':1-ar, 'pa':pa, 'x_center':0, 'y_center':0}

    lens_pert = {'Type':'pert', 'shear':shear, 'pa':pa_shear,\
                 'convergence':convergence, 'x_center':0, 'y_center':0}
    source_pos = {'x_center':xqso, 'y_center':yqso}

    lens = sim.OneObject(zl, mass_component=[lens_sie, lens_pert],\
               light_component=[])


    source = sim.OneObject(zs, light_component=[], positions=[source_pos])

    lenssystem = sim.LensSystem(lens, source)

    # the prefix of the GLAFIC output
    prefix = 'out'

    # search box: +/-10 * theta_E_max. should be good for almost all cases
    theta_max = 4 * np.pi * (sigma/3e5)** 2 * 206265

    x_range = (-10*theta_max, 10*theta_max)
    y_range = (-10*theta_max, 10*theta_max)

    minrange = np.min([x_range[1]-x_range[0], y_range[1]-y_range[0]])

    # pixel size: 0.01 theta_E_max. Sufficient to resolve the caustics
    pix_poi = 0.05 * theta_max

    imagedat = lenssystem.findimage(pix_poi=pix_poi, x_range=x_range, y_range=y_range)
    return imagedat

def lensing_crude_catalog(tbl, output=None):
    '''
    Run GLAFIC for all the rows in tbl adn summarize the info

    Parameters:
        tbl: astropy.table.Table
            The crude catalog generated by save_crude_catalog

        output: string
            The filename to save th lensed catalog.

    Return: tbl
        The lensed catalog
    '''

    # Quantities to save
    Nimg = []
    image_x_list = []
    image_y_list = []
    image_mu_list = []
    image_dt_list = []
    sep = []
    theta_E_list = []

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
        if zs > zl+0.01:
            image_x, image_y, image_mu, image_dt = \
                lens_one_object(qra, qdec, gra, gdec,\
                                zs, zl, sigma, ar, pa,\
                               convergence, shear, pa_shear)

        else:
            image_x, image_y, image_mu, image_dt = \
                    [np.array([-999.])] * 4

        image_x_list.append(fill_with_default(image_x,5))
        image_y_list.append(fill_with_default(image_y,5))
        image_mu_list.append(fill_with_default(image_mu,5))
        image_dt_list.append(fill_with_default(image_dt,5))

        Nimg.append(len(image_x))

        if zs > zl:
            theta_E_list.append(theta_E(sigma, zl, zs))
            xlist = image_x[image_x>-100]
            ylist = image_y[image_y>-100]
            sep.append(maxsep(xlist, ylist))

        else:
            theta_E_list.append(-999.)
            sep.append(-999.)

    tbl['image_x'] = np.array(image_x_list)
    tbl['image_y'] = np.array(image_y_list)
    tbl['image_mu'] = np.array(image_mu_list)
    tbl['image_dt'] = np.array(image_dt_list)
    tbl['Nimg'] = np.array(Nimg)
    tbl['maxsep'] = np.array(sep)
    tbl['theta_E'] = np.array(theta_E_list)

    tbl = tbl[tbl['Nimg']>1]

    if output is not None:
        tbl.write(output, overwrite=True)

    return tbl


