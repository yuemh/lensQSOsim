import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, vstack
from astropy.io import fits

import qsopop as qp
import utils
import lensing


def mockqso(zrange, outdir, abslim=-21, **kwargs):

    # the area is scaled to 1/16 of the sky
    allsky_deg2 = 4 * np.pi * (180/np.pi)**2
    cosmodc2_deg2 = 440.
    random_deg2 = allsky_deg2 / 16.

    # one realization = 1/50 sky
    area = allsky_deg2 / cosmodc2_deg2 * random_deg2 * 0.02
    # use a small area to quickly test the code
    area = 10 # for test

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

    # lens it
    tbl_lens = lensing.lensing_crude_catalog(crudecat)

    # save catalog
    if output is not None:
        tbl_lens.write(output, overwrite=True)

    return tbl_lens

def generate_mainsample(experiment, step):

    '''
    This is an example of generating a lensed quasar catalog.
    Please change parameters (mainly some paths) before running this function.

    I divide the code into two steps; step 1 is the generation of the mock quasars
    and step 2 is the lensing part. You can choose the step to run so that the function
    has more flexibility.
    '''

    abslim = -21 # the faint limit of M1450

    # define files and directories

    # the CosmoDC2 galaxy info file
    galinfofile = '/Users/minghao/Research/Projects/lensQSOsim/data/CosmoDC2/allmassivegalaxies_v1_short.fits'

    if experiment == 'highz':
        # the position to store the mock quasars
        mockqsodir = '/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/defaultbright/'

        # the redshift ranges to simulate; see qsopop.py for detailed information.
        # It is not very straightforward because SIMQSO takes observed magnitudes rather than absolute magnitudes as input
        zrangelist = ['z5', 'z6', 'z7']

        # the final output mock quasar catalog and the mock spectra
        qsoinfofile = mockqsodir + '/mockQSO_default_info_highz.fits'
        qsospecfile = mockqsodir + '/mockQSO_default_spec_highz.fits'

        # the position to store the realizations of the lensed quasars
        outputdir = '/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/test/highz'
        # the file to save the lensed quasar info
        outputfile = 'highzlensqso.fits'

        # Number of random realizations. One realization is 1/50 sky
        total_real = 10

    elif experiment == 'lowz':
        zrangelist = ['lowz', 'medz', 'z4']
        mockqsodir = '/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/test/'
        qsoinfofile = mockqsodir + '/mockQSO_default_info_lowz.fits'
        qsospecfile = mockqsodir + '/mockQSO_default_spec.fits'
        outputdir = '/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/test/lowz/'
        outputfile = 'lowzaddlensqso.fits'

        total_real = 50

    '''
    Custom paths definition end.
    No need to change the following code.
    '''

    # step 1: generate mock quasars
    if step == 1:
        os.system('mkdir -p %s'%mockqsodir)

        for zrange in zrangelist:
            mockqso(zrange, mockqsodir, abslim)

        # combine all the files

        combinedinfo = []
        combinedspec = []
        for zrange in zrangelist:
            infotbl = Table.read(mockqsodir+'/%sqso_M21.fits'%zrange)
            combinedinfo.append(infotbl)

            spechdu = fits.open(mockqsodir+'/%sqso_spectra.fits'%zrange)
            combinedspec.append(spechdu[0].data[infotbl['qidx']])

        combinedinfotbl = vstack(combinedinfo)
        combinedspecarr = np.concatenate(combinedspec)

        templatehdulist = spechdu.copy()
        templatehdulist[0].data = combinedspecarr
        templatehdulist[0].header['MAXIS2'] = len(combinedspecarr)
        combinedinfotbl['qidx'] = range(len(combinedinfotbl))

        templatehdulist.writeto(qsospecfile)
        combinedinfotbl.write(qsoinfofile)

    # step 2: lens the quasars
    if step == 2:
        os.system('mkdir -p %s'%outputdir)

        output = outputdir + '/' + outputfile

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
        # this is a temporary output
        tbl_allreals.write(outputdir+'/allr.fits', overwrite=True)

        # then, clean up column names
        tbl_allreals.remove_columns(['slopes', 'emLines', 'igmlos', 'synMag', 'synFlux',\
                                     'obsMag', 'obsMagErr', 'obsFlux', 'obsFluxErr', 'appMag'])
        tbl_allreals.rename_column('z', 'zqso')
        tbl_allreals.rename_column('redshift_true', 'zgal')
        tbl_allreals.rename_column('ra', 'gra')
        tbl_allreals.rename_column('dec', 'gdec')

        tbl_allreals.write(outputdir+'/'+outputfile)


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

def main():
#    generate_mainsample('highz', step=1)
#    generate_mainsample('highz', step=2)
#    addphot()
#    generate_simimage()

if __name__=='__main__':
    main()
