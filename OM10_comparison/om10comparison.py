import os, sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table, hstack, vstack
from astropy.io import fits, ascii
from astropy.cosmology import FlatLambdaCDM

defaultcosmo = FlatLambdaCDM(H0=70, Om0=0.3)


'''
Run this function first to generate an all-sky OM10 mock catalog.
'''
def select_OM10_allsky():
    '''
    Selecting an all-sky sample from the OM10 catalog.
    Requires om10- to be installed.
    '''
    # please change the path to your OM10 directory

    import om10
    db = om10.DB(catalog='/Users/minghao/Research/Projects/lensQSO/OM10/OM10/data/qso_mock.fits')
    db.select_random(maglim=60.0,area=41253,IQ=0.5)
    db.sample.write('./om10_allsky_sample.fits')


def plot_OM10_comparison_counts():
    # load data. Please change the data paths to your directories
    tbl_OM10 = Table.read('./om10_allsky_sample.fits')# this is OM10
    tbl_matchz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lensqso_allsky_allinfo.fits')# this is Yue+22

    # Now calculate the absolute magnitude for OM10 quasars
    # the MABS column seems to be wrong. So we re-calculate it.
    # Use Richards+06 K-correction. File attached.
    richards2006kcorr = ascii.read('richards2006kcorr.txt', format='cds')

    Miz2 = []

    cosmo_OM10 = FlatLambdaCDM(H0=72, Om0=0.26)
    for index in range(len(tbl_OM10)):
        redshift = tbl_OM10['ZSRC'][index]# quasar redshift
        # MAGI_IN is the intrinsic (un-magnified) apparent magnitude
        Miz2.append(tbl_OM10['MAGI_IN'][index] - cosmo_OM10.distmod(redshift).value - kcorr(redshift, richards2006kcorr))# this is Mi(z=2).

    tbl_OM10['Miz2'] = Miz2

    # next, select subsamples for each table
    # note that M1450 = Mi(z=2) + 1.486 (Ross+13, appendix B)
    tbl_OM10_select = tbl_OM10[tbl_OM10['Miz2']<-24-1.486]
    tbl_thiswork_select = tbl_matchz[(tbl_matchz['absMag']<-24)&(tbl_matchz['maxsep']>0.5)]


    plt.hist(tbl_thiswork_select['z'], bins=np.arange(0, 5.6, 0.2), label='This Work', alpha=0.4)
    plt.hist(tbl_OM10_select['ZSRC'], bins=np.arange(0, 5.6, 0.2), label='OM10', alpha=0.4)
    plt.xlabel('Quasar Redshift', fontsize=16)
    plt.ylabel(r'$N_{lens}$ (All Sky)', fontsize=16)
    plt.title(r'$M_{1450}<-24, \Delta \theta>$0.5"', fontsize=18)
    plt.legend(fontsize=14)
    plt.xlim([0., 3])
    plt.savefig('./OM10comparison_counts.pdf')
    plt.show()

def kcorr(z, tbl):
    zlist = tbl['z']
    kclist = tbl['KCorr']

    return np.interp(x=z, xp=zlist, fp=kclist)



plot_OM10_comparison_counts()
