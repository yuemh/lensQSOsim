import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, vstack
from astropy.io import fits
sys.path.append('/Users/minghao/Research/Projects/lensQSO/code2/dependencies/simqso')

from simqso import qsoSimulation,lumfun,sqmodels
dr9cosmo = FlatLambdaCDM(70,1-0.7,name='BOSSDR9')

import lensQSOsim.qsopop as qp
import lensQSOsim.glaficlensing as glafic
import lensQSOsim.utils as utils
import lensQSOsim.lensing as lensing


def mockqso(zrange, outdir, abslim=-21, **kwargs):

    # the area is scaled to 1/16 of the sky
    allsky_deg2 = 4 * np.pi * (180/np.pi)**2
    cosmodc2_deg2 = 440.
    random_deg2 = allsky_deg2 / 16.

    # one realization = 1/50 sky
    area = allsky_deg2 / cosmodc2_deg2 * random_deg2 * 0.02
#    area = 10 # for test

    # the name of the output file
    outname = zrange + 'qso'

    if "Ms" in kwargs.keys() and "beta" in kwargs.keys() and "alpha" in kwargs.keys():
        Ms = kwargs['Ms']
        beta = kwargs['beta']
        alpha = kwargs['alpha']
        outname += '_Ms%.1f_beta%.1f_alpha%.1f'%(Ms, beta, alpha)
        simParams = qp.make_simparams(zrange, area, outname, seed=12345, Mlim=abslim,\
                                      Ms=Ms, beta=beta, alpha=alpha)

    else:
        simParams = qp.make_simparams(zrange, area, outname, seed=12345, Mlim=abslim)

    qsoSimulation(simParams, saveSpectra=True, verbose=1, nproc=4)

    os.system('mv %s.fits %s'%(outname, outdir))
    os.system('mv %s_spectra.fits %s'%(outname, outdir))

    # selecting objects with desired M1450
    tbl = Table.read(outdir+'/'+outname+'.fits')
    tbl['qidx'] = range(len(tbl))

    mask = (tbl['absMag']<abslim)

    tbl_bright = tbl[mask]
    tbl_bright.write(outdir+'/'+outname+'_M%d.fits'%(-abslim), overwrite=True)

def lens_one_realization(qsoinfofile, galinfofile,\
                         output=None, seed=12345):
    np.random.seed(seed)

    # assign random coordinates to quasars

    qsoinfo = Table.read(qsoinfofile)
    qra, qdec = utils.random_positions(ra_range=[45, 90],\
                                       dec_range=[-90,0],\
                                       Ncoord=len(qsoinfo))

    qsoinfo['qra'] = qra
    qsoinfo['qdec'] = qdec

    # load galaxy info
    galinfo = Table.read(galinfofile)

    # find random matched between mock quasars and mock galaxies
    crudecat = lensing.save_crude_catalog(qsoinfo, galinfo, seed=seed)

    crudecat.write('./crudecattest.fits', overwrite=True)
    # lens it
    tbl_lens = lensing.lensing_crude_catalog(crudecat)

    # save catalog
    if output is not None:
        tbl_lens.write(output, overwrite=True)

    return tbl_lens

def generate_mainsample(experiment, step):
    abslim = -21
    galinfofile = '/Users/minghao/Research/Projects/lensQSOsim/data/CosmoDC2/allmassivegalaxies_v1_short.fits'

    if experiment == 'highz':
        mockqsodir = '/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/default/'
        qsoinfofile = mockqsodir + '/mockQSO_default_info_highz.fits'
        outputdir = '/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/highz'

        total_real = 1000

    elif experiment == 'lowz':
        mockqsodir = '/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/default/'
        qsoinfofile = mockqsodir + '/mockQSO_default_info_lowz.fits'
        outputdir = '/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lowz'

        total_real = 50

    # step 1: generate mock quasars
    if step == 1:
        os.system('mkdir -p %s'%mockqsodir)

        if experiment == 'highz':
            mockqso('z4', mockqsodir, abslim)
            mockqso('z5', mockqsodir, abslim)
            mockqso('z6', mockqsodir, abslim)
            mockqso('z7', mockqsodir, abslim)

        elif experiment == 'lowz':
            mockqso('lowz', mockqsodir, abslim)
            mockqso('medz', mockqsodir, abslim)

    # step 2: lens the quasars
    if step == 2:
        os.system('mkdir -p %s'%outputdir)

        output = outputdir + '/allr.fits'

        allreals = []
        for index in range(0, total_real, 1):
            if os.path.exists(outputdir + '/r%d.fits'%index):
                allreals.append(Table.read(outputdir + '/r%d.fits'%index))
                continue

            tbl_lens = lens_one_realization(qsoinfofile, galinfofile, seed=index)

            #print(tbl_lens)
            if len(tbl_lens)>0:
                tbl_lens['realization'] = index
                tbl_lens.write(outputdir + '/r%d.fits'%index)

                allreals.append(tbl_lens)

        tbl_allreals = vstack(allreals)
        tbl_allreals.write(output, overwrite=True)


def z6paramtest(step=2):

    galinfofile = '/Users/minghao/Research/Projects/lensQSOsim/data/CosmoDC2/allmassivegalaxies_v1_short.fits'

    Ms = -24.90
    alpha = -1.23
    for beta in [-3.4, -3.0, -2.6]:

        mockqsodir = '/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/z6paramtest/var_alpha'
        qsoinfofile = mockqsodir + '/z6slopetestqso_Ms%.1f_beta%.1f_alpha%.1f_M21.fits'%(Ms, beta, alpha)
        outputdir = '/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/z6paramtest' + \
                '/z6slopetestqso_Ms%.1f_beta%.1f_alpha%.1f'%(Ms, beta, alpha)

        total_real = 200

        if step==1:
            mockqso('z6slopetest', mockqsodir, -21, Ms=Ms, beta=beta, alpha=alpha)

        elif step==2:
            os.system('mkdir -p %s'%outputdir)

            output = outputdir + '/allr.fits'

            allreals = []
            for index in range(0, total_real, 1):
                if os.path.exists(outputdir + '/r%d.fits'%index):
                    allreals.append(Table.read(outputdir + '/r%d.fits'%index))
                    continue

                tbl_lens = lens_one_realization(qsoinfofile, galinfofile, seed=index)

                #print(tbl_lens)
                if len(tbl_lens)>0:
                    tbl_lens['realization'] = index
                    tbl_lens.write(outputdir + '/r%d.fits'%index)

                    allreals.append(tbl_lens)

            tbl_allreals = vstack(allreals)
            tbl_allreals.write(output, overwrite=True)


#    ztest = 3
#    Ntest = 100000
#    zlist = np.ones(Ntest) * ztest
#    ztesttbl = Table({'z': zlist, 'qidx': range(Ntest)})
#    ztesttbl.write(mockqsodir+'/z%dfixed_sis.fits'%ztest)


#    Ms = -24.90
#    beta = -2.73
#    for alpha in np.arange(-2, -1, 0.1):
#        mockqso('z6slopetest', mockqsodir, -21, Ms=Ms, beta=beta, alpha=alpha)

#    Ms = -24.90
#    alpha = -1.23
#    for beta in np.arange(-3.6, -2.6, 0.1):
#        mockqso('z6slopetest', mockqsodir, -21, Ms=Ms, beta=beta, alpha=alpha)

#    mockqso('z7a2', mockqsodir, -21)

#    Ms = -24.90
#    alpha = -1.23
#    for beta in np.arange(-3.6, -2.6, 0.1):
#        continue


def main():
#    generate_mainsample('highz', step=2)
    z6paramtest(step=2)

if __name__=='__main__':
    main()
