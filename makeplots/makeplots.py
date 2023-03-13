import os, sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table, hstack, vstack
from astropy.io import fits, ascii

from astropy.cosmology import FlatLambdaCDM
from scipy.special import gamma as gammafunc
from scipy.integrate import quad
from matplotlib.ticker import FormatStrFormatter

sys.path.append('/Users/minghao/Research/Projects/lensQSOsim/code/lensQSOsim')
from lensQSOsim import sed

defaultcosmo = FlatLambdaCDM(H0=70, Om0=0.3)

import matplotlib
import mylib.utils.zscale as zs

def plot_observable_counts():

    tbl_lowz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lowzlensqso_allinfo.fits')
    tbl_highz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/highzlensqso_allinfo.fits')

    print(np.max(tbl_highz['absMag']))


    # define subsets
    PS1_subset_lowz = get_discoverable(tbl_lowz, magindex=[3, 4], maglim=[23.1, 22.3], seplim=[1.11, 1.07])
    LSST_subset_lowz = get_discoverable(tbl_lowz, magindex=[3, 4], maglim=[26.8, 26.1], seplim=[0.71, 0.69])
    Legacy_subset_lowz = get_discoverable(tbl_lowz, magindex=[2, 4], maglim=[23.5, 22.5], seplim=[1.18, 1.11])
    DES_subset_lowz = get_discoverable(tbl_lowz, magindex=[3, 4], maglim=[23.8, 23.1], seplim=[0.88, 0.83])
    HSC_subset_lowz = get_discoverable(tbl_lowz, magindex=[3, 4], maglim=[26.4, 25.5], seplim=[0.56, 0.63])

    PS1_subset_highz = get_discoverable(tbl_highz, magindex=[3, 4], maglim=[23.1, 22.3], seplim=[1.11, 1.07])
    LSST_subset_highz = get_discoverable(tbl_highz, magindex=[3, 4], maglim=[26.8, 26.1], seplim=[0.71, 0.69])
    Legacy_subset_highz = get_discoverable(tbl_highz, magindex=[2, 4], maglim=[23.5, 22.5], seplim=[1.18, 1.11])
    DES_subset_highz = get_discoverable(tbl_highz, magindex=[3, 4], maglim=[23.8, 23.1], seplim=[0.88, 0.83])
    HSC_subset_highz = get_discoverable(tbl_highz, magindex=[3, 4], maglim=[26.4, 25.5], seplim=[0.56, 0.63])

    # report the numbers
    # note: high-z part corresponds to 20 all-sky

    print('LSST', 0.05*20000./40000.*len(LSST_subset_highz) + 20000./40000.*len(LSST_subset_lowz))
    print('PS1', 0.05*30000./40000.*len(PS1_subset_highz) + 30000./40000.*len(PS1_subset_lowz))
    print('legacy', 0.05*14000./40000.*len(Legacy_subset_highz) + 14000./40000.*len(Legacy_subset_lowz))
    print('DES', 0.05*5000./40000.*len(DES_subset_highz) + 5000./40000.*len(DES_subset_lowz))
    print('HSC', 0.05*1400./40000.*len(HSC_subset_highz) + 1400./40000.*len(HSC_subset_lowz))

    # plot

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=[12, 5])

    axes[0].hist(LSST_subset_lowz['z'], bins=np.arange(0, 5.2, 0.2), weights=[0.5]*len(LSST_subset_lowz), label='LSST')
    axes[0].hist(PS1_subset_lowz['z'], bins=np.arange(0, 5.2, 0.2), weights=[0.75]*len(PS1_subset_lowz), label='PS1')

    axes[1].hist(LSST_subset_highz['z'], bins=np.arange(5, 8, 0.2), weights=[0.05*0.5]*len(LSST_subset_highz), label='LSST')
    axes[1].hist(PS1_subset_highz['z'], bins=np.arange(5, 8, 0.2), weights=[0.05*0.75]*len(PS1_subset_highz), label='PS1')


    axes[0].set_xlim([0, 5])
    axes[1].set_xlim([5, 8])
    axes[1].legend()


    axes[0].set_xlabel(r'$z_{qso}$', fontsize=14)
    axes[1].set_xlabel(r'$z_{qso}$', fontsize=14)
    axes[0].set_ylabel(r'$N_{lens, discoverable}$', fontsize=14)
    axes[1].set_ylabel(r'$N_{lens, discoverable}$', fontsize=14)

    axes[0].set_title('Low Redshift (20 Realizations)', fontsize=14)
    axes[1].set_title('High Redshift (1,000 Realizations)', fontsize=14)

    plt.tight_layout()
    plt.savefig('./plots/estimatedcounts.pdf')


def plot_OM10_comparison():
    # load data
    tbl_OM10 = Table.read('./om10_allsky_sample.fits')

    tbl_matchz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lensqso_allsky_allinfo.fits')

    # select M1450 and sep
    # first, assign M1450 to OM10
    richards2006kcorr = ascii.read('./richards2006kcorr.txt', format='cds')

    Miz2 = []

    cosmo_OM10 = FlatLambdaCDM(H0=72, Om0=0.26)
    for index in range(len(tbl_OM10)):
        redshift = tbl_OM10['ZSRC'][index]
        Miz2.append(tbl_OM10['MAGI_IN'][index] - cosmo_OM10.distmod(redshift).value - kcorr(redshift, richards2006kcorr))

    tbl_OM10['Miz2'] = Miz2

    # next, select subsamples for each table
    tbl_OM10_select = tbl_OM10[tbl_OM10['Miz2']<-24-1.486]
    tbl_thiswork_select = tbl_matchz[(tbl_matchz['absMag']<-24)&(tbl_matchz['maxsep']>0.5)]


    plt.hist(tbl_thiswork_select['z'], bins=np.arange(0, 5.6, 0.2), label='This Work', alpha=0.4)
    plt.hist(tbl_OM10_select['ZSRC'], bins=np.arange(0, 5.6, 0.2), label='OM10', alpha=0.4)
    plt.xlabel('Quasar Redshift', fontsize=16)
    plt.ylabel(r'$N_{lens}$ (All Sky)', fontsize=16)
    plt.title(r'$M_{1450}<-24, \Delta \theta>$0.5"', fontsize=18)
    plt.legend(fontsize=14)
    plt.xlim([0., 3])
    plt.savefig('./plots/OM10comparison_counts.pdf')
    plt.show()

    return 0
    # plot QLFs

    mslist = np.arange(-30, -22, 0.1)
    zlist = [1, 2, 3]
    clist = ['black', 'purple', 'b']
    shiftlist = [+2, 0, -1]

    lines_to_legend_tw = []
    lines_to_legend_OM10 = []

    zlabels = ['z=%.1f'%z for z in zlist]
    worklabels = ['This Work', 'OM10']

    for index, z in enumerate(zlist):

        lf_om10 = np.log10(lf_OM10(mslist, z)) +shiftlist[index]
        lf_tw = np.log10(lf_dr9(mslist, z)) +shiftlist[index]

        lines_to_legend_tw.append(plt.plot(mslist, lf_tw, color=clist[index])[0])
        lines_to_legend_OM10.append(plt.plot(mslist, lf_om10, '--', color=clist[index])[0])


    plt.xlabel(r'$M_{1450}$', fontsize=14)
    plt.ylabel(r'$\log~\Phi~\mathdefault{[Mpc^{-3}Mag^{-1}]}$', fontsize=14)
    plt.xlim(-22,-30)

    plt.text(-23, -3.6, 'z=1 (+2 dex)', color=clist[0], fontsize=14)
    plt.text(-23, -5.4, 'z=2', color=clist[1], fontsize=14)
    plt.text(-23, -7.9, 'z=3 (-1 dex)', color=clist[2], fontsize=14)

    #plt.ylim([1e-8, 1e-6])
    plt.title('QLFs', fontsize=16)
    plt.legend(lines_to_legend_tw, zlabels, fontsize=14)
    plt.legend([lines_to_legend_tw[0], lines_to_legend_OM10[0]], worklabels, fontsize=14)
    plt.savefig('./plots/OM10comparison_QLF.pdf')
    plt.show()

    # plot VDFs

    meanlogsig_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/dc2meanlogsig.dat')
    logvdf_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/dc2logvdf.dat')
    logvdferr_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/dc2logvdferr.dat')


    sigmalist = np.arange(100, 400, 1)
    logvdf_om10 = np.log10(vdf_om10(sigmalist))
    plt.plot(np.log10(sigmalist), logvdf_om10, 'k--', label='OM10')

    colorlist = ['m', 'c', 'orange', 'r']

    logvdf, logvdferr, meanlogsig = np.loadtxt('../galpop/localvdf_CosmoDC2.dat')
    mask = (meanlogsig>2)
    plt.plot(meanlogsig[mask], logvdf[mask], '-', label='CosmoDC2, z<0.1', color='k')

    plt.fill_between(x=meanlogsig[mask], y1=logvdf[mask]-logvdferr[mask], y2=logvdf[mask]+logvdferr[mask],\
                   color='k', alpha=0.2)


    for index in range(len(logvdf_dc2)):

        zmin = 0.3 * (1+index)
        zmax = zmin + 0.3

        print(zmin, zmax)

        zindex = index

        meanlogsig = meanlogsig_dc2[zindex]
        logvdf = logvdf_dc2[zindex]
        logvdferr = logvdferr_dc2[zindex]

        plt.plot(meanlogsig, logvdf, '-', label='CosmoDC2, %.1f<z<%.1f'%(zmin, zmax), color=colorlist[index])
        plt.fill_between(x=meanlogsig, y1=logvdf-logvdferr, y2=logvdf+logvdferr,\
                   color=colorlist[index], alpha=0.2)

    plt.legend(fontsize=14)
    plt.xlim([1.95, 2.65])
    plt.xlabel(r'$\log \sigma$ [km/s]', fontsize=14)
    plt.ylabel(r'$\log \Phi~\mathdefault{[Mpc^3~dex^{-1}]}$', fontsize=14)
    plt.title('VDFs', fontsize=16)
    plt.savefig('./plots/OM10comparison_VDF.pdf')
    plt.show()

def plot_separation():
    tbl_matchz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lensqso_allsky_allinfo.fits')

    fig, ax = plt.subplots(figsize=[6,4])

    ax.hist(tbl_matchz['maxsep'], density=True, bins=np.arange(0, 5, 0.1), cumulative=True, histtype='step', edgecolor='k')
    ax.plot([0.15, 0.15], [-0.2, 1.2], 'k--')


    ax.set_xlabel(r'$\Delta \theta$ (arcsec)', fontsize=14)
    ax.set_ylabel(r'$P(<\Delta \theta)$', fontsize=14)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 4])
    plt.tight_layout()
    plt.savefig('separation.pdf')


def plot_zd():

    tbl_matchz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lensqso_allsky_allinfo.fits')

    tbl_z2 = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/fixed_z/z2fixed/allr.fits')
    tbl_z5 = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/fixed_z/z5fixed/allr.fits')

    tbl_allqso = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/default/mockQSO_default_info.fits')
    tbl_allqso = tbl_allqso[tbl_allqso['absMag']<-20+0.89]

    print(np.median(tbl_matchz['redshift_true']))

    fig, ax = plt.subplots(figsize=[6,4])

    weight_mockq = [1./50./440./(len(tbl_allqso)/(41253./16.))/0.1] * len(tbl_matchz)
    weight_fixz2 = [1./50./440./(100000/(41253./16.))/0.1] * len(tbl_z2)
    weight_fixz5 = [1./50./440./(100000/(41253./16.))/0.1] * len(tbl_z5)

    ax.hist(tbl_matchz['redshift_true'], facecolor='k', edgecolor='None', bins=np.arange(0,3.1,0.1), label='All', weights=weight_mockq, alpha=0.5)
    ax.hist(tbl_z2['redshift_true'], histtype='step', edgecolor='b', bins=np.arange(0,3.1,0.1), label=r'$z_s=2$', weights=weight_fixz2, lw=2)
    ax.hist(tbl_z5['redshift_true'], histtype='step', edgecolor='r', bins=np.arange(0,3.1,0.1), label=r'$z_s=5$', weights=weight_fixz5, lw=2)


    ax.set_xlabel('Deflector Redshift', fontsize=12)
    ax.legend()


    ax.set_ylabel(r'$dP/dz$', fontsize=12)

#    plt.tight_layout()
    plt.savefig('./plots/zd.pdf')
    plt.show()

def plot_example(index, qsospechdu, galinfotbl, selectidx=None, comment='', sample='lowz'):
    tbl_lowz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/%s/allrphot.fits'%sample)
    tbl_lowz['allidx'] = range(len(tbl_lowz))

    tbl = add_quantities(tbl_lowz)

    if selectidx is None:

        tbl = tbl[(tbl['Qmag'][:,4]<26)&(tbl['Gmag'][:,4]<26)&(tbl['maxsep']>0.4)\
                  &(tbl['secondmag'][:,4]-tbl['firstmag'][:,4]<1)]

    else:
        tbl = tbl[tbl['allidx']==selectidx]
        index = 0

    ### get the spectra of the quasar and the galaxy

#    index = 1
    print(tbl['allidx', 'z', 'redshift_true', 'maxsep'][index])
    print(tbl['image_mu'][index])

    # quasar
    image_mu = tbl['image_mu'][index]
    mask = image_mu > -900
    image_mu = image_mu[mask]
    image_mu = np.abs(image_mu)
    image_mu_sorted = np.sort(image_mu)


    qidx = tbl['qidxnew'][index]
#    qsospechdu = fits.open('/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/default/mockQSO_default_spec.fits')[0]
    qsospechdr = qsospechdu.header
    qsospecdat = qsospechdu.data

    qsowave = np.exp(qsospechdr['CRVAL1']\
                  + np.arange(qsospechdr['NAXIS1']) * qsospechdr['CD1_1'])
    qsospec = qsospecdat[qidx]*1e-17*np.sum(image_mu)

    # galaxy
    gidx = tbl['galaxy_id'][index]
#    galinfotbl = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/CosmoDC2/allmassivegalaxies_v1.fits')

#    galinfotbl = tbl

    sedkeys = []
    for col in galinfotbl.colnames:
        if col.startswith('sed') and\
           not ('disk' in col or 'bulge' in col or 'ext' in col):
            sedkeys.append(col)

    sedkeys = sed.sort_sed_keys(sedkeys)

    redshift = tbl['redshift_true'][index]
    lumdist = defaultcosmo.luminosity_distance(redshift).to('cm').value / 1e15

    # tophat spec
    row = galinfotbl[galinfotbl['galaxy_id']==gidx]
    sedlist = [row[key][0] for key in sedkeys]
    thspec = sed.tophat_cosmoDC2_spec(sedlist, sedkeys, redshift)

    # template spec
    tpspec = sed.template_cosmoDC2_spec(sedlist, redshift)

    # re-grid

    galwave = qsowave
    mask_th = (galwave<np.max(thspec.wavelength)+50)

    galwave_th = galwave[mask_th]
    galwave_tp = galwave[~mask_th]

    galflux_th = np.interp(galwave_th, thspec.wavelength, thspec.value)
    galflux_tp = np.interp(galwave_tp, tpspec.wavelength, tpspec.value)

    galspec = np.concatenate([galflux_th, galflux_tp])


    # photometries
    mags = tbl['apermag'][index]

    wavelengths = [3680.04, 4782.26, 6217.82, 7532.28, 8685.06, 9730.07,\
                   12483.00, 16313.00, 22010.00, 33526.00, 46028.00]

    fnu_mags = 3631 * 10 ** (-0.4 * mags) * 1e-23
    flambda_mags = fnu_mags * 3e18 / np.array(wavelengths)**2

    simimagedir = '/Users/minghao/Research/Projects/lensQSOsim/data/simimage/default/%s/'%sample

    rimage = fits.open(simimagedir+'/LSSTr/tmp/nsy%d.fits'%selectidx)[0].data
    iimage = fits.open(simimagedir+'/LSSTi/tmp/nsy%d.fits'%selectidx)[0].data
    zimage = fits.open(simimagedir+'/LSSTz/tmp/nsy%d.fits'%selectidx)[0].data

#    vmin, vmax = zs.zscale(zimage)

    vmin = None
#    vmax = np.max(iimage)

    rrgb = raw_to_rgb(rimage, vmin=vmin, vmax=np.max(rimage)/1.2)
    irgb = raw_to_rgb(iimage, vmin=vmin, vmax=np.max(iimage)/1.2)
    zrgb = raw_to_rgb(zimage, vmin=vmin, vmax=np.max(zimage)/1.2)

    rgb = np.moveaxis(np.array([zrgb, irgb, rrgb]), 0, -1)

    plt.close('all')

    ### plot everything
    plt.close('all')
    fig, axes = plt.subplots(figsize=[10,4], gridspec_kw={'width_ratios': [1,1.5]}, ncols=2)

    axes[1].plot(qsowave, qsospec, 'b-' ,alpha=0.5, label='Quasar (Magnified)')
    axes[1].plot(galwave, galspec, 'r-', alpha=0.5, label='Galaxy')
    axes[1].plot(galwave, qsospec+galspec, 'k-', lw=2, alpha=0.5, label='All')
    axes[1].plot(wavelengths, flambda_mags, 'ko', label='Synthetic Magnitudes')

    axes[1].set_xlabel('Wavelength [Angstrom]', fontsize=14)
    axes[1].set_ylabel(r'Flux $\mathdefault{[erg~s^{-1}~cm^{-2}~\AA^{-1}]}$', fontsize=14)
    axes[1].legend(fontsize=14)

    axes[1].set_title(r'$z_{qso}=%.2f, z_d=%.2f$'%(tbl['z'][index], tbl['redshift_true'][index]), fontsize=14)

    axes[1].set_xscale('log')
    axes[1].set_xticks([5e3, 1e4, 5e4])
    axes[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axes[1].tick_params(axis='x', labelsize=12)
    axes[1].tick_params(axis='y', labelsize=12)

    axes[0].axis('off')

    axes[0].imshow(rgb, origin='lower')

    axes[0].set_title('LSST $riz$', fontsize=16)

    axes[0].plot([5,10], [5,5], 'w-')
    axes[0].text(6, 6, '1"', fontsize=14, color='w')

    axes[0].text(3, 42, comment, fontsize=14, color='w')

    plt.tight_layout()
    plt.savefig('./examples/%s/example_%d.pdf'%(sample,selectidx))
    plt.close('all')


def plot_quad():
    allqso = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/default/mockQSO_default_info.fits')
    area = 4834.678442744569

    tbl = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lensqso_allsky_allinfo.fits')

    LSST_resolve = tbl[tbl['maxsep']>0.666*0.7]
    LSST_quads = LSST_resolve[LSST_resolve['Nimg']>3]

#    LSST_quads_bright = LSST_quads[LSST_quads['thirdmag'][:,4]<26.1]
#    print(len(LSST_quads_bright))

    # get the number of quads
    tbl_lowz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lowzlensqso_allinfo.fits')
    tbl_highz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/highzlensqso_allinfo.fits')

    PS1_subset_lowz = get_discoverable(tbl_lowz, magindex=[3, 4], maglim=[23.1, 22.3], seplim=[1.11, 1.07])
    LSST_subset_lowz = get_discoverable(tbl_lowz, magindex=[3, 4], maglim=[26.8, 26.1], seplim=[0.71, 0.69])
    Legacy_subset_lowz = get_discoverable(tbl_lowz, magindex=[2, 4], maglim=[23.5, 22.5], seplim=[1.18, 1.11])
    DES_subset_lowz = get_discoverable(tbl_lowz, magindex=[3, 4], maglim=[23.8, 23.1], seplim=[0.88, 0.83])
    HSC_subset_lowz = get_discoverable(tbl_lowz, magindex=[3, 4], maglim=[26.4, 25.5], seplim=[0.56, 0.63])

    PS1_subset_lowz_quad = PS1_subset_lowz[PS1_subset_lowz['Nimg']>3]
    LSST_subset_lowz_quad = LSST_subset_lowz[LSST_subset_lowz['Nimg']>3]
    Legacy_subset_lowz_quad = Legacy_subset_lowz[Legacy_subset_lowz['Nimg']>3]
    DES_subset_lowz_quad = DES_subset_lowz[DES_subset_lowz['Nimg']>3]
    HSC_subset_lowz_quad = HSC_subset_lowz[HSC_subset_lowz['Nimg']>3]

    PS1_subset_highz = get_discoverable(tbl_highz, magindex=[3, 4], maglim=[23.1, 22.3], seplim=[1.11, 1.07])
    LSST_subset_highz = get_discoverable(tbl_highz, magindex=[3, 4], maglim=[26.8, 26.1], seplim=[0.71, 0.69])
    Legacy_subset_highz = get_discoverable(tbl_highz, magindex=[2, 4], maglim=[23.5, 22.5], seplim=[1.18, 1.11])
    DES_subset_highz = get_discoverable(tbl_highz, magindex=[3, 4], maglim=[23.8, 23.1], seplim=[0.88, 0.83])
    HSC_subset_highz = get_discoverable(tbl_highz, magindex=[3, 4], maglim=[26.4, 25.5], seplim=[0.56, 0.63])

    PS1_subset_highz_quad = PS1_subset_highz[PS1_subset_highz['Nimg']>3]
    LSST_subset_highz_quad = LSST_subset_highz[LSST_subset_highz['Nimg']>3]
    Legacy_subset_highz_quad = Legacy_subset_highz[Legacy_subset_highz['Nimg']>3]
    DES_subset_highz_quad = DES_subset_highz[DES_subset_highz['Nimg']>3]
    HSC_subset_highz_quad = HSC_subset_highz[HSC_subset_highz['Nimg']>3]

    print('LSST', 0.05*20000./40000.*len(LSST_subset_highz_quad) + 20000./40000.*len(LSST_subset_lowz_quad))
    print('PS1', 0.05*30000./40000.*len(PS1_subset_highz_quad) + 30000./40000.*len(PS1_subset_lowz_quad))
    print('legacy', 0.05*14000./40000.*len(Legacy_subset_highz_quad) + 14000./40000.*len(Legacy_subset_lowz_quad))
    print('DES', 0.05*5000./40000.*len(DES_subset_highz_quad) + 5000./40000.*len(DES_subset_lowz_quad))
    print('HSC', 0.05*1400./40000.*len(HSC_subset_highz_quad) + 1400./40000.*len(HSC_subset_lowz_quad))

    # non-lense
    PS1_nonlens = allqso[(allqso['synMag'][:,3]<23.1)|(allqso['synMag'][:,4]<22.3)]
    LSST_nonlens = allqso[(allqso['synMag'][:,3]<26.8)|(allqso['synMag'][:,4]<26.1)]
    Legacy_nonlens = allqso[(allqso['synMag'][:,3]<23.5)|(allqso['synMag'][:,4]<22.5)]
    DES_nonlens = allqso[(allqso['synMag'][:,3]<23.8)|(allqso['synMag'][:,4]<23.1)]
    HSC_nonlens = allqso[(allqso['synMag'][:,3]<26.4)|(allqso['synMag'][:,4]<25.5)]

    print('LSST', 20000./area*len(LSST_nonlens))
    print('PS1', 30000./area*len(PS1_nonlens))
    print('legacy', 14000./area*len(Legacy_nonlens))
    print('DES', 5000./area*len(DES_nonlens))
    print('HSC', 1400./area*len(HSC_nonlens))



    fig, axes = plt.subplots(ncols=2, figsize=[10,4])

    axes[0].hist(allqso['synMag'][:,3], bins=np.arange(18, 29, 1), histtype='step', label='unlensed',\
                weights=[20000/area]*len(allqso), edgecolor='gray', alpha=0.5)

    axes[0].hist(LSST_resolve['criteriamag'][:,3], bins=np.arange(18, 29, 1), histtype='step', edgecolor='k',\
                label='lensed', weights=[0.5]*len(LSST_resolve))
    axes[0].hist(LSST_quads['criteriamag'][:,3], bins=np.arange(18, 29, 1), facecolor='gray', label='quad',\
                weights=[0.5]*len(LSST_quads))
    axes[0].set_yscale('log')
    axes[0].set_xlabel('LSST $i$', fontsize=12)
    axes[0].set_ylabel('$N_{lens}$', fontsize=12)

    axes[1].set_xlabel('$e_{galaxy}$', fontsize=12)
    axes[1].set_ylabel('$dP/de$', fontsize=12)

#    axes[0].legend(loc='upper left', bbox_to_anchor=[0,1])
    axes[0].legend()
#    axes[0].set_ylim([1e-2, 1e8])

    gauss_mean = 1-0.7
    gauss_sigma = 0.16

    arlist = np.arange(0.0, 1.01, 0.01)

    problist = 1/gauss_sigma/np.sqrt(2*np.pi) * np.exp(-(arlist-gauss_mean)**2/2/gauss_sigma**2)

    tbl_cosmodc2 = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/CosmoDC2/allmassivegalaxies_v1_short_correctq.fits')
    elist = 1 - tbl_cosmodc2['axisratio']

    mean_cosmodc2 = np.mean(elist)

    # get the mean for om10
    func1 = lambda x: 1/gauss_sigma/np.sqrt(2*np.pi) * np.exp(-(x-gauss_mean)**2/2/gauss_sigma**2)*x
    func0 = lambda x: 1/gauss_sigma/np.sqrt(2*np.pi) * np.exp(-(x-gauss_mean)**2/2/gauss_sigma**2)

    mean_om10 = quad(func1, 0, 1)[0] / quad(func0, 0, 1)[0]

    axes[1].hist(elist, histtype='step', bins = np.arange(0,1.01,0.01), density=True, label='CosmoDC2', edgecolor='k')
    axes[1].plot(arlist, problist, 'k--', label='OM10')

    axes[1].set_xlim([0, 1])

    axes[1].legend()

    plt.tight_layout()
    plt.savefig('./plots/quad.pdf')
    plt.show()
    print(mean_om10, mean_cosmodc2)


def plot_td():

    tbl = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lensqso_allsky_allinfo.fits')

    LSST_subset = tbl[tbl['maxsep']>0.666*0.71]

    fig, axes = plt.subplots(ncols=4, figsize=[13,4], sharey=True)
    plt.subplots_adjust(left=0.06, right=0.99, bottom=0.15, hspace=0.01, wspace=0.03)

    # double
    LSST_double = LSST_subset[LSST_subset['Nimg']<3]
    mu1d, dt1d = get_dmag_dt(LSST_double, imgidx=0)
    mu2d, dt2d = get_dmag_dt(LSST_double, imgidx=1)

    mag1d = -2.5 * np.log10(mu1d)
    mag2d = -2.5 * np.log10(mu2d)

    axes[0].plot(mag1d-mag2d, dt2d-dt1d, 'k.', ms=0.5)
    axes[0].set_yscale('log')

    # quad
    LSST_quad = LSST_subset[LSST_subset['Nimg']>3]

    mu1, dt1 = get_dmag_dt(LSST_quad, imgidx=0)
    mu2, dt2 = get_dmag_dt(LSST_quad, imgidx=1)
    mu3, dt3 = get_dmag_dt(LSST_quad, imgidx=2)
    mu4, dt4 = get_dmag_dt(LSST_quad, imgidx=3)

    mag1 = -2.5 * np.log10(mu1)
    mag2 = -2.5 * np.log10(mu2)
    mag3 = -2.5 * np.log10(mu3)
    mag4 = -2.5 * np.log10(mu4)

    axes[1].plot(mag1-mag2, dt2-dt1, 'k.', ms=3)
    axes[1].set_yscale('log')
    axes[2].plot(mag1-mag3, dt3-dt1, 'k.', ms=3)
    axes[2].set_yscale('log')
    axes[3].plot(mag1-mag4, dt4-dt1, 'k.', ms=3)
    axes[3].set_yscale('log')

    axes[0].set_title('Double, x=2', fontsize=12)
    axes[1].set_title('Quad, x=2', fontsize=12)
    axes[2].set_title('Quad, x=3', fontsize=12)
    axes[3].set_title('Quad, x=4', fontsize=12)

    fig.supxlabel('$\Delta m=m_1-m_x$ [mags]', fontsize=12)
    axes[0].set_ylabel('$\Delta t=t_x-t_1$ [days]', fontsize=12)

    plt.tight_layout()
    plt.savefig('./plots/dmdt.pdf')
    plt.show()

    mag_2bright_d = get_mu_nth(LSST_double, nth=2)
    mag_3bright_q = get_mu_nth(LSST_quad, nth=3)

    print(np.mean(mag_2bright_d))
    print(np.mean(mag_3bright_q))

def plot_deflector_flux():

    tbl = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lensqso_allsky_allinfo.fits')

    LSST_subset = tbl[((tbl['criteriamag'][:,3]<26.8)&(tbl['maxsep']>0.666*0.71))\
                    |((tbl['criteriamag'][:,4]<26.1)&(tbl['maxsep']>0.666*0.69))]

    fig, axes = plt.subplots(ncols=2, figsize=[10,4])

    # plot the flux
    axes[0].hist(LSST_subset['galmags'][:,3], histtype='step', edgecolor='k', bins=np.arange(14, 28, 0.5), label='LSST discoverable', weights=[0.5]*len(LSST_subset))
    axes[0].set_xlabel(r'$m_{i,galaxy}$', fontsize=14)
    axes[0].set_ylabel(r'$N_{lens}$ (LSST Area)', fontsize=14)
    axes[0].legend(loc='upper left')
    axes[0].plot([26.8, 26.8], [0,270], 'k--')
    axes[0].set_ylim([0, 270])

    print(tbl.colnames)

    axes[1].set_xlabel(r'$m_{i,galaxy}-m_{i,quasar}$', fontsize=14)
    axes[1].hist(LSST_subset['Gmag'][:,3]-LSST_subset['Qmag'][:,3], histtype='step', edgecolor='k', bins=np.arange(-10, 10, 0.5), label='LSST discoverable', weights=[0.5]*len(LSST_subset))
    axes[1].plot([0, 0], [0,270], 'k--')
    axes[1].set_ylim([0, 220])

    axes[0].text(24, 250, 'LSST limit')
    axes[1].text(-5.5, 190, 'Galaxy\nDominated')
    axes[1].text(2.5, 190, 'Quasar\nDominated')

    plt.tight_layout()
    plt.savefig('./plots/deflectorflux.pdf')
    plt.show()


'''
Below are some utilities
'''

def add_quantities(tbl):
    Qmag = []
    Gmag = []
    apermag = []
    firstmag = []
    secondmag = []
    thirdmag = []
    fourthmag = []
    criteriamag = []

    for index in range(len(tbl)):
        image_mu = tbl['image_mu'][index]
        mask = image_mu > -900
        image_mu = image_mu[mask]
        image_mu = np.abs(image_mu)
        image_mu_sorted = np.sort(image_mu)

        qsomag_raw = tbl['qsomag'][index]
        qsomag_magnified = qsomag_raw - 2.5 * np.log10(np.sum(image_mu_sorted))
        Qmag.append(qsomag_magnified)

        qsomag_first = qsomag_raw - 2.5 * np.log10(image_mu_sorted[-1])
        firstmag.append(qsomag_first)

        qsomag_second = qsomag_raw - 2.5 * np.log10(image_mu_sorted[-2])
        secondmag.append(qsomag_second)

        if len(image_mu_sorted)>3:
            qsomag_third = qsomag_raw - 2.5 * np.log10(image_mu_sorted[-3])
            thirdmag.append(qsomag_third)

            qsomag_fourth = qsomag_raw - 2.5 * np.log10(image_mu_sorted[-4])
            fourthmag.append(qsomag_fourth)

            criteriamag.append(qsomag_third)

        else:
            thirdmag.append([-999.]*len(qsomag_second))
            fourthmag.append([-999.]*len(qsomag_second))
            criteriamag.append(qsomag_second)

        galmag = tbl['galmags'][index]
        Gmag.append(galmag)

        allmag = -2.5 * np.log10(10**(-0.4*galmag)+10**(-0.4*qsomag_magnified))
        apermag.append(allmag)

    tbl['Qmag'] = np.array(Qmag)
    tbl['Gmag'] = np.array(Gmag)
    tbl['apermag'] = np.array(apermag)
    tbl['firstmag'] = np.array(firstmag)
    tbl['secondmag'] = np.array(secondmag)
    tbl['thirdmag'] = np.array(thirdmag)
    tbl['fourthmag'] = np.array(fourthmag)
    tbl['criteriamag'] = np.array(criteriamag)

    return tbl


def kcorr(z, tbl):
    zlist = tbl['z']
    kclist = tbl['KCorr']

    return np.interp(x=z, xp=zlist, fp=kclist)


def lf_OM10(M, z):

    h = 0.72
    a = 2.98
    b = 4.05
    zs = 1.6

    fz = np.exp(a*z) * (1+np.exp(b*zs)) / (np.sqrt(np.exp(b*z)) + np.sqrt(np.exp(b*zs)))**2

    # this is M1450
    Ms = -20.9 + 5 * np.log10(h) - 2.5 * np.log10(fz) + 0.89

    Phis = 5.34e-6 * h**3

    if z>3:
        beta, alpha = -1.45, -2.58
    else:
        beta, alpha = -1.45, -3.31

    lf = Phis / (10**(0.4*(alpha+1)*(M-Ms))+10**(0.4*(beta+1)*(M-Ms)))

    return lf

def lf_dr9_ple(M, z):
    alpha, beta = -1.16, -3.37
    k1, k2 = 1.241, -0.249
    Phis = 10** -5.96
    Ms = -22.85 - 2.5 * (k1 * z + k2 * z**2) + 1.486

    lf = Phis / (10**(0.4*(alpha+1)*(M-Ms))+10**(0.4*(beta+1)*(M-Ms)))

    return lf

def lf_dr9_lede(M, z):
    alpha, beta = -1.29, -3.51
    c1, c2 = -0.689, -0.809

    logPhis = -5.93 + c1 * (z-2.2)
    Ms = -26.57 + c2 * (z-2.2) + 1.486

    Phis = 10** logPhis

    lf = Phis / (10**(0.4*(alpha+1)*(M-Ms))+10**(0.4*(beta+1)*(M-Ms)))

    return lf

def lf_dr9(M, z):
    if z<2.2:
        return lf_dr9_ple(M, z)
    elif z>2.2:
        return lf_dr9_lede(M, z)

def vdf_om10(sigma):
    h = 0.72

    Phis = 8e-3 * h **3 # Mpc^-3
    ss = 161
    alpha, beta = 2.32, 2.67

    vdf = Phis * (sigma/ss)**alpha \
        * np.exp(-(sigma/ss)**beta) \
        * beta / gammafunc(alpha/beta) * np.log(10)

    # Mpc^-3 dex^-1
    return vdf

def schechter1(v, phis, vs, alpha):
    term1 = (v/vs)**alpha
    term2 = np.exp(-(v/vs)**4)
    term3 = 4 * (v/vs)**4 * np.log(10)
    
    return phis * term1 * term2 * term3

def vdf_wyithe(sigma):
    phis = 0.27e-2
    alpha = -1.19*4
    vs = 220
    #vs = 161
    
    return schechter1(sigma, phis, vs, alpha)


def get_discoverable(tbl, magkey='criteriamag', magindex=[3,4], maglim=[28, 28], seplim=[0.7, 0.7]):
    mask = [False] * len(tbl)

    for ifilt in range(len(magindex)):
        mi = magindex[ifilt]
        mlim = maglim[ifilt]
        sep = seplim[ifilt]
        mask = mask|((tbl[magkey][:,mi]<mlim)&(tbl['maxsep']>sep*0.666))

    return tbl[mask]

def raw_to_rgb(image, vmin=None, vmax=None):
    if vmin is None:
        vmin = np.min(image)

    if vmax is None:
        vmax = np.max(image)

    rgb = (image-np.min(image))/(vmax-np.min(image))
    rgb[rgb<0]=0
    rgb[rgb>1]=1
    
    return rgb

def testcosmodc2():

    meanlogsig_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/cosmodc2vdfs/allmeanlogsig.dat')
    logvdf_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/cosmodc2vdfs/alllogvdf.dat')
    logvdferr_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/cosmodc2vdfs/alllogvdferr.dat')

    sigmalist = np.arange(100, 400, 1)
    logvdf_om10 = np.log10(vdf_om10(sigmalist))
    plt.plot(np.log10(sigmalist), logvdf_om10, 'k--', label='OM10')

    colorlist = ['m', 'c', 'orange', 'r']

    logvdf, logvdferr, meanlogsig = np.loadtxt('../galpop/localvdf_CosmoDC2.dat')
    plt.plot(meanlogsig, logvdf, '-', label='CosmoDC2, z<0.1', color='k')

    plt.fill_between(x=meanlogsig, y1=logvdf-logvdferr, y2=logvdf+logvdferr,\
                   color='k', alpha=0.2)


    for index in range(10):
        zindex = index
        zmin = zindex * 0.1 + 0.05
        zmax = zmin + 0.1

        print(zmin, zmax)
        meanlogsig = meanlogsig_dc2[index]
        logvdf = logvdf_dc2[index]
        logvdferr = logvdferr_dc2[index]

        plt.plot(meanlogsig, logvdf, '-', label='CosmoDC2, %.1f<z<%.1f'%(zmin, zmax))
        plt.fill_between(x=meanlogsig, y1=logvdf-logvdferr, y2=logvdf+logvdferr, alpha=0.2)

#    logsigmalist, logphi, logphiep, logphiem = np.loadtxt('./localvdf.dat')
#    plt.errorbar(x=logsigmalist, y=logphi, yerr=[logphiem, logphiep], fmt='o',\
#                    color='k', label='Hasan+2019', mfc='None')


    plt.legend(fontsize=14)
    plt.xlim([1.95, 2.65])
    plt.xlabel(r'$\log \sigma$ [km/s]', fontsize=14)
    plt.ylabel(r'$\log \Phi~\mathdefault{[Mpc^3~dex^{-1}]}$', fontsize=14)
    plt.title('VDFs', fontsize=16)
#    plt.savefig('./plots/OM10comparison_VDF.pdf')
    plt.show()


def get_dmag_dt(tbl, imgidx):
    mulist = []
    dtlist = []

    for index in range(len(tbl)):
        mu = tbl['image_mu'][index]
        dt = tbl['image_dt'][index]

        mask = (mu>-900)
        mu = mu[mask]
        dt = dt[mask]

        mu = np.abs(mu)
        sortindex = np.argsort(dt)

        mu_sort = mu[sortindex]
        dt_sort = dt[sortindex]

        mulist.append(mu_sort[imgidx])
        dtlist.append(dt_sort[imgidx])

    return np.array(mulist),\
            np.array(dtlist)


def get_mu_nth(tbl, nth=1):
    mulist = []

    for index in range(len(tbl)):
        mu = tbl['image_mu'][index]

        mask = (mu>-900)
        mu = mu[mask]

        mu = np.abs(mu)
        sortindex = np.argsort(mu)

        mu_sort = mu[sortindex][::-1]
#        print(mu_sort)
#        break

        mulist.append(mu_sort[nth-1])

    return np.array(mulist)

def save_tbl_allinfo():
    tbl_lowz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lowz/allrphot.fits')
    tbl_highz = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/highz/allrphot.fits')

    tbl_lowz = add_quantities(tbl_lowz)
    tbl_lowz.write('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lowzlensqso_allinfo.fits', overwrite=True)

    tbl_highz = add_quantities(tbl_highz)
    tbl_highz.write('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/highzlensqso_allinfo.fits', overwrite=True)

    tbl = vstack([tbl_lowz, tbl_highz[(tbl_highz['realization']<50)]])
    tbl.write('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/lensqso_allsky_allinfo.fits', overwrite=True)


# old, unused plotting functions
def plot_detectable_frac():
    tbl_b36_lens = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/z6paramtest/z6slopetestqso_Ms-24.9_beta-2.6_alpha-1.2/allrphot.fits')
    tbl_b36_all = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/z6paramtest/var_alpha/z6slopetestqso_Ms-24.9_beta-3.4_alpha-1.2_M21.fits')
    tbl_b36_lens = add_quantities(tbl_b36_lens)

#    print(tbl_b36_all.meta)

    tbl_b32_lens = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/z6paramtest/z6slopetestqso_Ms-24.9_beta-3.0_alpha-1.2/allrphot.fits')
    tbl_b32_all = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/z6paramtest/var_alpha/z6slopetestqso_Ms-24.9_beta-3.0_alpha-1.2_M21.fits')
    tbl_b32_lens = add_quantities(tbl_b32_lens)

    tbl_b28_lens = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/z6paramtest/z6slopetestqso_Ms-24.9_beta-2.6_alpha-1.2/allrphot.fits')
    tbl_b28_all = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/z6paramtest/var_alpha/z6slopetestqso_Ms-24.9_beta-2.6_alpha-1.2_M21.fits')
    tbl_b28_lens = add_quantities(tbl_b28_lens)

    maglims = np.arange(20.2, 23, 0.2)
    n36lens_discover = []
    n36lens_detect = []
    n36all = []

    n32lens_discover = []
    n32lens_detect = []
    n32all = []

    n28lens_discover = []
    n28lens_detect = []
    n28all = []

    for m in maglims:
        n36lens_discover.append(len(tbl_b36_lens[(tbl_b36_lens['criteriamag'][:,4]<m)&(tbl_b36_lens['maxsep']>0.7*0.666)]))
        n36lens_detect.append(len(tbl_b36_lens[(tbl_b36_lens['firstmag'][:,4]<m)]))
        n36all.append(len(tbl_b36_all[tbl_b36_all['synMag'][:,4]<m]))

        n32lens_discover.append(len(tbl_b32_lens[(tbl_b32_lens['criteriamag'][:,4]<m)&(tbl_b32_lens['maxsep']>0.7*0.666)]))
        n32lens_detect.append(len(tbl_b32_lens[(tbl_b32_lens['firstmag'][:,4]<m)]))
        n32all.append(len(tbl_b32_all[tbl_b32_all['synMag'][:,4]<m]))

        n28lens_discover.append(len(tbl_b28_lens[(tbl_b28_lens['criteriamag'][:,4]<m)&(tbl_b28_lens['maxsep']>0.7*0.666)]))
        n28lens_detect.append(len(tbl_b28_lens[(tbl_b28_lens['firstmag'][:,4]<m)]))
        n28all.append(len(tbl_b28_all[tbl_b28_all['synMag'][:,4]<m]))


    tau = (len(tbl_b36_lens)+len(tbl_b32_lens)+len(tbl_b28_lens))/20/41253/(len(tbl_b36_all)+len(tbl_b32_all)+len(tbl_b28_all))\
            * (41253/16.)
    print(tau)

    np.savetxt('./detectfrac.txt', [n36lens_detect, n36all, n32lens_detect, n32all, n28lens_detect, n28all])
    return 0

    fig, ax = plt.subplots()

    ax.errorbar(maglims, np.array(n28lens_detect)/20/(np.array(n28all)/24173.392213722847*41253),\
                 yerr=np.sqrt(n28lens_detect)/20/(np.array(n28all)/24173.392213722847*41253),\
                 fmt='o-', color='black', label=r'$\beta=-2.6$')
    ax.errorbar(maglims, np.array(n32lens_detect)/20/(np.array(n32all)/24173.392213722847*41253),\
                 yerr=np.sqrt(n32lens_detect)/20/(np.array(n32all)/24173.392213722847*41253),\
                 fmt='o-', color='darkorange', label=r'$\beta=-3.0$')
    ax.errorbar(maglims, np.array(n36lens_detect)/20/(np.array(n36all)/24173.392213722847*41253),\
                 yerr=np.sqrt(n36lens_detect)/20/(np.array(n36all)/24173.392213722847*41253),\
                 fmt='o-', color='red', label=r'$\beta=-3.4$')

#    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    ax.plot([21.8, 21.8], [-0.005, 0.05], 'k--')
    ax.text(21.2, 0.017, r'$M_{1450}=M_*$', fontsize=12)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.legend(fontsize=12)
    ax.set_xlabel(r'$m_{z,lim}$', fontsize=14)
    ax.set_ylabel(r'$N_{lens}~/~N_{all}$', fontsize=14)
    ax.set_title(r'Detectable', fontsize=14)

    ax.set_ylim([0, 0.04])

    plt.tight_layout()
    plt.savefig('./plots/detect_frac.pdf')

    plt.show()


    fig, ax = plt.subplots()

    ax.errorbar(maglims, np.array(n28lens_discover)/20/(np.array(n28all)/24173.392213722847*41253),\
                 yerr=np.sqrt(n28lens_discover)/20/(np.array(n28all)/24173.392213722847*41253),\
                 fmt='o-', color='black', label=r'$\beta=-2.8$')
    ax.errorbar(maglims, np.array(n32lens_discover)/20/(np.array(n32all)/24173.392213722847*41253),\
                 yerr=np.sqrt(n32lens_discover)/20/(np.array(n32all)/24173.392213722847*41253),\
                 fmt='o-', color='darkorange', label=r'$\beta=-3.2$')
    ax.errorbar(maglims, np.array(n36lens_discover)/20/(np.array(n36all)/24173.392213722847*41253),\
                 yerr=np.sqrt(n36lens_discover)/20/(np.array(n36all)/24173.392213722847*41253),\
                 fmt='o-', color='red', label=r'$\beta=-3.6$')

#    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    ax.plot([21.8, 21.8], [-0.005, 0.02], 'k--')
    ax.text(21.2, 0.003, r'$M_{1450}=M_*$', fontsize=12)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.legend(fontsize=12)
    ax.set_xlabel(r'$m_{z,lim}$', fontsize=14)
    ax.set_ylabel(r'$N_{lens}~/~N_{all}$', fontsize=14)
    ax.set_title(r'Discoverable', fontsize=14)

    ax.set_ylim([0, 0.007])

    plt.tight_layout()
    plt.savefig('./plots/discover_frac.pdf')
    plt.show()


def plot_vdf_wyithe():

    sigmalist = np.arange(100, 400, 1)
#    logvdftest = np.log10(vdf_om10(sigmalist))
    logvdf_wyithe =  np.log10(vdf_wyithe(sigmalist))

#    logsigmalist, logphi, logphiep, logphiem = np.loadtxt('./localvdf.dat')

#    fig, ax = plt.subplots(figsize=[12,8])
#    plt.plot(np.log10(sigmalist), logvdftest, label='OM10')
#    plt.plot(np.log10(sigmalist), logvdfwyithe, label='Wyithe+02')
#    plt.plot(logv[mask], Flogv[mask], 'g-', label='Comerford+02, SIS')
#    plt.plot(logv[~mask], Flogv[~mask], 'g--', label='Comerford+02, NFW')
#    plt.fill_between(x=logsigmalist, y1=logphi-logphiem, y2=logphi+logphiep,\
#                        color='gray', alpha=0.2, label='local')


    meanlogsig_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/cosmodc2vdfs/allmeanlogsig.dat')
    logvdf_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/cosmodc2vdfs/alllogvdf.dat')
    logvdferr_dc2 = np.loadtxt('/Users/minghao/Research/Projects/lensQSOsim/code/galpop/cosmodc2vdfs/alllogvdferr.dat')

    sigmalist = np.arange(100, 400, 1)
    plt.plot(np.log10(sigmalist), logvdf_wyithe, 'k--', label='Wyithe+2002')


    logvdf, logvdferr, meanlogsig = np.loadtxt('../galpop/localvdf_CosmoDC2.dat')
    mask = (meanlogsig>2)
    plt.plot(meanlogsig[mask], logvdf[mask], '-', label='CosmoDC2, z<0.1', color='k')

    plt.fill_between(x=meanlogsig[mask], y1=logvdf[mask]-logvdferr[mask], y2=logvdf[mask]+logvdferr[mask],\
                   color='k', alpha=0.2)


    colorlist = ['m', 'c', 'orange', 'r']

    for index in range(4):
        zindex = 4 + 5 * index
        zmin = zindex * 0.1 + 0.05
        zmax = zmin + 0.1
        meanlogsig = meanlogsig_dc2[zindex]
        logvdf = logvdf_dc2[zindex]
        logvdferr = logvdferr_dc2[zindex]

        plt.plot(meanlogsig, logvdf, '-', label='CosmoDC2, %.2f<z<%.2f'%(zmin, zmax), color=colorlist[index])
        plt.fill_between(x=meanlogsig, y1=logvdf-logvdferr, y2=logvdf+logvdferr,\
                   color=colorlist[index], alpha=0.2)


    plt.xlabel(r'$\log \sigma$ [km/s]', fontsize=14)
    plt.ylabel(r'$\log \Phi~\mathdefault{[Mpc^3~dex^{-1}]}$', fontsize=14)
    plt.legend(fontsize=12)
    plt.xlim([1.95, 2.65])
    plt.ylim([-6, -1])
    plt.title('VDFs', fontsize=14)
    plt.tight_layout()

    plt.savefig('./plots/wyithevdf.pdf')
    plt.show()


def plot_yang17():
    tbl = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/highz/allrphot_yang17.fits')

    LSST_subset = tbl[((tbl['criteriamag'][:,5]<24.9)&(tbl['maxsep']>0.666*0.67))\
                    |((tbl['criteriamag'][:,4]<26.1)&(tbl['maxsep']>0.666*0.67))]


    LSST_Qmag = LSST_subset[LSST_subset['qselect']]
    LSST_QGmag = LSST_subset[LSST_subset['qgselect']]
    LSST_QGmagmorph = LSST_QGmag[LSST_QGmag['dmag']<0.5]

    plt.hist(LSST_subset['z'], bins=np.arange(5.0, 6.0, 0.1), label='All LSST discoverable',\
            edgecolor='k', histtype='step')
    plt.hist(LSST_Qmag['z'], bins=np.arange(5.0, 6.0, 0.1), label='Q-Only Color',\
            edgecolor='None', facecolor='silver')
    plt.hist(LSST_QGmag['z'], bins=np.arange(5.0, 6.0, 0.1), label='Q+D Color',\
            edgecolor='None', facecolor='steelblue')
    plt.hist(LSST_QGmagmorph['z'], bins=np.arange(5.0, 6.0, 0.1), label='Q+D Color and Morphology',\
            edgecolor='None', facecolor='darkorange')

    plt.legend()
    plt.xlabel('Redshift', fontsize=12)
    plt.ylabel('Number', fontsize=12)

    plt.title('Candidate Selection Completeness', fontsize=12)
    plt.savefig('./plots/completeness.pdf')
    plt.show()


def plot_ugmag():
    tbl = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/mocklens/default/highz/allrphot_yang17.fits')

    print(len(tbl))

    LSST_subset = tbl[((tbl['criteriamag'][:,5]<24.9)&(tbl['maxsep']>0.666*0.67))\
                    |((tbl['criteriamag'][:,4]<26.1)&(tbl['maxsep']>0.666*0.67))]

#    plt.hist(LSST_subset['Gmag'][:,0], edgecolor='purple', histtype='step',
#             bins=np.arange(15, 35, 0.1))
    plt.hist(LSST_subset['Gmag'][:,1], edgecolor='k', histtype='step',
             bins=np.arange(15, 35, 0.2), label='LSST Discoverable')

    plt.plot([27.4,27.4], [0,250], 'k--')
    plt.text(28, 100, 'LSST limit', fontsize=12)
    plt.ylim([0, 130])
    plt.xlabel('Deflector g-band magnitude', fontsize=12)
    plt.ylabel('Number', fontsize=12)
    plt.legend()
    plt.savefig('./plots/dropout.pdf')
    plt.show()


# the main function

def main():
#    plot_observable_counts()
#    plot_OM10_comparison()
#    plot_separation()
    plot_zd()
#    plot_quad()
#    plot_td()
#    plot_deflector_flux()

    '''
    qsospechdu = fits.open('/Users/minghao/Research/Projects/lensQSOsim/data/mockQSO/default/mockQSO_default_spec.fits')[0]
    galinfotbl = Table.read('/Users/minghao/Research/Projects/lensQSOsim/data/CosmoDC2/allmassivegalaxies_v1.fits')


    plot_example(0, qsospechdu=qsospechdu, galinfotbl=galinfotbl, selectidx=251, comment='High-z Quad', sample='highz')
    plot_example(0, qsospechdu=qsospechdu, galinfotbl=galinfotbl, selectidx=228, comment='High-z Double\nQuasar-Dominated', sample='highz')
    plot_example(0, qsospechdu=qsospechdu, galinfotbl=galinfotbl, selectidx=442, comment='High-z Double\nDeflector-Dominated', sample='highz')
    plot_example(0, qsospechdu=qsospechdu, galinfotbl=galinfotbl, selectidx=385, comment='High-z Double\nMarginally Discoverable', sample='highz')

    plot_example(0, qsospechdu=qsospechdu, galinfotbl=galinfotbl, selectidx=7, comment='Low-z Quad', sample='lowz')
    plot_example(0, qsospechdu=qsospechdu, galinfotbl=galinfotbl, selectidx=268, comment='Low-z Double\nQuasar-Dominated', sample='lowz')
    plot_example(0, qsospechdu=qsospechdu, galinfotbl=galinfotbl, selectidx=131, comment='Low-z Double\nDeflector-Dominated', sample='lowz')
    plot_example(0, qsospechdu=qsospechdu, galinfotbl=galinfotbl, selectidx=309, comment='Low-z Double\nMarginally Discoverable', sample='lowz')
    '''
#    testcosmodc2()
#    save_tbl_allinfo()


#    plot_detectable_frac()
#    plot_vdf_wyithe()
#    plot_yang17()
#    plot_ugmag()

#    plot_detectable_frac()


if __name__=='__main__':
    main()
