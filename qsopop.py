#!/usr/bin/env python

import os, sys
import numpy as np
from astropy.cosmology import FlatLambdaCDM
defaultcosmo = FlatLambdaCDM(H0=70, Om0=0.3)

import matplotlib.pyplot as plt


from setpaths import simqsodir
sys.path.append(simqsodir)
from simqso import qsoSimulation,lumfun,sqmodels


'''
The following two QLFs are from Ross+2013. 
lowzQLF: SDSS DR9 PLE, for z<2.2
medzQLF: SDSS DR9 LEDE, for 2.2<z<3.5
'''

def lowzQLF(which=1):
    # PLE model in Ross+2013, for 0<z<2.2
    if which==1:
        row = -1.16,-3.37,-22.85,1.241,-0.249,-5.96
    alpha,beta,MStar_i_z0,k1,k2,logPhiStar = row
    MStar1450_z0 = MStar_i_z0 + 1.486
    MStar = lumfun.PolyEvolParam([-2.5*k2,-2.5*k1,MStar1450_z0])
    return lumfun.DoublePowerLawLF(logPhiStar,MStar,alpha,beta,
                                   cosmo=defaultcosmo)

def medzQLF():
    # LEDE model in Ross+2013, for 2.2<z<4
    c1,c2 = -0.689, -0.809
    logPhiStar_z2_2 = -5.93
    MStar_i_z2_2 = -26.57
    MStar1450_z0 = MStar_i_z2_2 + 1.486 # --> M1450
    MStar = lumfun.PolyEvolParam([c2,MStar1450_z0],z0=2.2)
    logPhiStar = lumfun.PolyEvolParam([c1,logPhiStar_z2_2],z0=2.2)
    alpha = -1.29
    beta = -3.51
    return lumfun.DoublePowerLawLF(logPhiStar,MStar,alpha,beta,
                                   cosmo=defaultcosmo)

'''
The following QLF is based on Matsuoka+2018 (HSC z~6 QLF),
but the slopes will be replaced by the input.
Used to test the impact of bright-end slope
'''

def M18QLFvarslope(Ms, beta, alpha):
    # z~6 QLF in Matsuoka+2018, for 5.5<z<6.5, Standard (Table 5)
    c1,c2 = -0.7, 0
    logPhiStar_z6 = np.log10(10.9/1e9)
    MStar1450 = Ms # --> M1450
    MStar = lumfun.PolyEvolParam([c2,MStar1450], z0=6)
    logPhiStar = lumfun.PolyEvolParam([c1,logPhiStar_z6], z0=6)
    alpha = alpha
    beta = beta
    return lumfun.DoublePowerLawLF(logPhiStar,MStar,alpha,beta,
                                   cosmo=defaultcosmo)


'''
The following two QLFs are from Wang+2019, for z>6.5.
W19QLFb3: re-fit with beta ~ -3
W19QLFb4: re-fit with beta = -4
'''

def W19QLFb3():
    # z>6.7
    logPhiStar_z67 = -8.433509799321072
    logPhiStar = lumfun.PolyEvolParam([-0.78, logPhiStar_z67], z0=6.7)
    MStar = -24.9
    beta = -2.51
    alpha = -1.23
    return lumfun.DoublePowerLawLF(logPhiStar,MStar,alpha,beta,
                                   cosmo=defaultcosmo)


def W19QLForig():
    logPhiStar_z6 = np.log10(3.17e-9)
    logPhiStar = lumfun.PolyEvolParam([-0.78, logPhiStar_z6], z0=6.7)

    MStar = -25.2
    alpha, beta = -1.90, -2.54
    return lumfun.DoublePowerLawLF(logPhiStar,MStar,alpha,beta,
                                   cosmo=defaultcosmo)


'''
The following QLFs are the combined z>3.5 QLF with linear interpolation between redshifts.
z=3.5: DR9
z=3.9:Akiyama+18
z=5.0: Kim+20
z=6.0: Matsuoka+18
z=6.7: Wang+19
'''

def highzQLF_b3():
    zlist = [3.5, 3.9, 5, 6.0, 6.7]
    Mslist = [-26.1357, -25.36, -25.78, -24.9, -24.9]
    logPhilist = [-6.8257, -6.57511836, -7.36, -7.9625735, -8.433509799321072]
    alphalist = [-3.51, -3.11, -3.44, -2.73, -2.5114194345960246]
    betalist = [-1.29, -1.3, -1.21, -1.23, -1.23]

    logPhiStar = lumfun.LinearInterpEvolParam([zlist, logPhilist])
    MStar = lumfun.LinearInterpEvolParam([zlist, Mslist])
    alpha = lumfun.LinearInterpEvolParam([zlist, alphalist])
    beta = lumfun.LinearInterpEvolParam([zlist, betalist])

    return lumfun.DoublePowerLawLF(logPhiStar,MStar,alpha,beta,
                                   cosmo=defaultcosmo)

def make_simparams(zrange, area, outname, Mlim=-21, seed=12345, Ms=-25, beta=-2.5, alpha=-1.3):
    # determine the redshift range
    if zrange=='lowz':
        zmin = 0.1
        zmax = 2.2
        lumfun = lowzQLF()
        maglim = 24 + (Mlim+21)

    elif zrange == 'medz':
        zmin = 2.2
        zmax = 3.5
        lumfun = medzQLF()
        maglim = 26 + (Mlim+21)

    elif zrange == 'z4':
        zmin = 3.5
        zmax = 4.5
        lumfun = highzQLF_b3()
        maglim = 25.2 + (Mlim+21)

    elif zrange == 'z5':
        zmin = 4.5
        zmax = 5.5
        lumfun = highzQLF_b3()
        maglim = 25.7 + (Mlim+21)

    elif zrange == 'z6':
        zmin = 5.5
        zmax = 6.7
        lumfun = highzQLF_b3()
        maglim = 26.5 + (Mlim+21)

    elif zrange == 'z7':
        zmin = 6.7
        zmax = 7.5
        lumfun = W19QLFb3()
        maglim = 27 + (Mlim+21)

    ###########################
    #####   Test slopes   #####
    ###########################

    elif zrange == 'z6slopetest':
        zmin = 5.5
        zmax = 6.5

        lumfun = M18QLFvarslope(Ms=Ms, beta=beta, alpha=alpha)
        maglim = 26.5 + (Mlim+21)

    else:
        raise KeyError('Unknown redshift range key %s'%zrange)

    simParams = {
      # filename for simulation output (".fits" is appended)
      'FileName':outname,
    #  'GridFileName':'gridtest',
      # wavelength range of simulated spectra in Angstroms
      'waveRange':(3000.,60000),
      # constant spectral dispersion per pixel as R = d(log(lam))
      'SpecDispersion':500,
      # dispersion scale is logarithmic [only option for now]
      'DispersionScale':'logarithmic',
      # set the cosmology, any appropriate instance from astropy.cosmology allowed
      'Cosmology':defaultcosmo,
      # setting a global random seed allows the simulation to be repeatable
      'RandomSeed':seed,
      # Define the "grid" of points in (M,z) space for the simulation
      # In this case the grid is a distribution of points sampled from the    
      # Ross et al. 2013 QLF determined from BOSS DR9.
      'GridParams':{
        # Define the grid as coming from a luminosity function
        'GridType':'LuminosityFunction', 
        # Specify the functional form of the LF, using a double power-law
        # with evolutionary parameters from Ross et al. 2013
        'QLFmodel': lumfun,
        # simulate a 10k deg2 survey
        #'QLFargs':{'skyArea':1e4},       
        'QLFargs':{'skyArea':area}, # but only 10 deg2 will run a lot faster
        # set bright and faint flux limits
        'mRange':(14, maglim),
        # and redshift range
        'zRange':(zmin, zmax),
        # flux range defined in r-band
          'ObsBand':'SDSS-z',
        # rest luminosities defined at 1450A
        'RestBand':1450.,
      },
      # Define the model for generating quasar emission spectra
      'QuasarModelParams':{
        # underlying continuum, only option is power law
        'ContinuumParams':{
          # power law slopes have a gaussian distribution
          'ContinuumModel':'BrokenPowerLaw',
          # the continuum consists of a series of broken power laws with
          # independent slope distributions 
          # the format is [ (meanSlope, stdSlope), breakWavelength, ... ]
          'PowerLawSlopes':[(-1.5,0.3),1100,(-0.5,0.3),
                        5700,(-0.37,0.3),9730,(-1.7,0.3),22300,(-1.03,0.3)],
        },
        # the emission line model
        'EmissionLineParams':{
          # the emission line profile distribution comes from the BOSS DR9 model
          # allowing for the Baldwin Effect
          'EmissionLineModel':'VariedEmissionLineGrid',
          # these are rescalings of the equivalent widths in the model, determined
          # emprically by matching colors with BOSS quasars. I.e., the Lyman alpha
          # EW is 10% greater than the nominal value in the model.
          'scaleEWs':{'LyAb':1.1,'LyAn':1.1,
                      'CIVb':0.75,'CIVn':0.75,
                      'CIII]b':0.8,'CIII]n':0.8,
                      'MgIIb':0.8,'MgIIn':0.8},
        },
        # the Fe emission template from Vestergaard & Wilkes 2001
        'IronEmissionParams':{
          # rescalings of sections of the template, empirically determined fitting
          # of composite BOSS quasar spectra
          'FeScalings':[(0,1540,0.5),(1540,1680,2.0),(1680,1868,1.6),
                        (1868,2140,1.0),(2140,3500,1.0)],
        },
      },
      # define the model for transmission spectra through the HI forest, based on
      # Monte Carlo realizations of absorption systems
      'ForestParams':{
        # filename to save the forest transmission spectra
    #    'FileName':'boss_dr9qlf_forest',
        # name of the model for the distribution of absorbers 
        'ForestModel':sqmodels.forestModels['McGreer+2013'], 
        # redshift range over which to sample absorbers
        'zRange':(0.0, 8.0),
        # the number of independent sightlines to generate
        # WP11 suggest a minimum of 2000 to properly sample the scatter
        'NumLinesOfSight':4000,
        #'NumLinesOfSight':200, # however, 200 will run a lot faster for testing
        # the minimum spectral dispersion to use when generating the transmission
        # spectra; R=30000 => 10 km/s is a good value to capture the weak systems
        'Rmin':30000.,
      },
      # define the photometric systems for the survey, namely, the bandpasses
      # for calculating synthetic photometry from the spectra, and an error model
      # for producing realistic fluxes and errors
      'PhotoMapParams':{
        # list the systems individually, the output 'synMag' and 'obsMag' fields
        # will have a final dimension equal to the total number of bandpasses,
        # in the order listed here. I.e., for this simulation the shape is
        # synMag[...,9], with synMag[...,0:5] representing SDSS 'ugriz' and
        # synMag[...,5:9] representing UKIDSS 'YJHK'
        'PhotoSystems':[
          # SDSS ugriz with an error model for single-epoch imaging from the
          # Legacy footprint. A subset of the bandpasses could be specified using,
          # e.g., ('SDSS','Legacy','griz'),
          ('SDSS','Legacy'),
#          ('LSST','Wide', 'y'),
          # UKIDSS YJHK with an error model for the Large Area Survey
#          ('UKIRT','UKIDSS_LAS', 'JHK'),
#          ('WISE','AllWISE'),
        ]
      },
      'maxFeatureIter': 1
    }

    return simParams



