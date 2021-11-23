import os, sys
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import astropy.constants as const

from scipy.integrate import quad, dblquad

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

'''
The analytical VDF
'''
def vdf_analytical(sigma, z):
    a = -0.15
    b = 2.35
    sigma0 = 172.2 * (1+z)**0.18

    phis_z0 = 5.86e-3 * np.e / np.log(10)

    phis = phis_z0 * (1+z) ** -1.18

    return phis * (sigma/sigma0) ** a * np.exp(-(sigma/sigma0)**b) * np.log(10)


'''
Functions for tau_m
'''

def theta_E_func(sigma, zl, zs, cosmo):

    Ds = cosmo.angular_diameter_distance(zs)
    Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)

    return 4 * np.pi * (sigma/3e5)**2 * Dls / Ds


def inverse_theta_E_func(theta_E, zl, zs, cosmo):
    Ds = cosmo.angular_diameter_distance(zs)
    Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)

    sc2 = theta_E / 4 / np.pi * Ds / Dls

    return np.sqrt(sc2) * 3e5

def tau_integral(sigma, z, zs, vdf, cosmo):

    if z>zs:
        return 0

    vdfterm = vdf(sigma, z) / sigma / np.log(10) * u.Mpc**-3

    dV = cosmo.differential_comoving_volume(z)

    theta_E = theta_E_func(sigma, z, zs, cosmo)
    area  = np.pi * theta_E **2  * u.sr

    return vdfterm * dV * area


def sep_integral(theta_E, z, zs, vdf, cosmo):
    if z>zs:
        return 0

    sigma = inverse_theta_E_func(theta_E, z, zs, cosmo)

#    print(sigma)
    vdfterm = vdf(sigma, z) / sigma / np.log(10) * u.Mpc**-3

    dV = cosmo.differential_comoving_volume(z)
    Ds = cosmo.angular_diameter_distance(zs)
    Dls = cosmo.angular_diameter_distance_z1z2(z, zs)

    theta_E = theta_E_func(sigma, z, zs, cosmo)
    area  = np.pi * theta_E **2  * u.sr
    additional_factor = 8 * np.pi * sigma / (3e5)**2 * Dls / Ds

    return vdfterm * dV * area / additional_factor


def taum(zs, vdf, cosmo=cosmo):
    paras = [zs, vdf, cosmo]
    result = dblquad(tau_integral, 0, zs, 0, 1000, args=paras)
    return result[0]


def taumdiff(zd, zs, vdf, cosmo=cosmo):
    paras = (zd, zs, vdf, cosmo)
    result = quad(tau_integral, 0, 1000, args=paras)

    return result[0]


def sep_distribution_diff(theta_E, zs, vdf, cosmo=cosmo):
    theta_E_rad = theta_E / 206265
    paras = (zs,  vdf, cosmo)

    intfunc = lambda z: sep_integral(theta_E_rad, z, zs, vdf, cosmo)

    result = quad(intfunc, 0, zs)

    return result[0] / 206265

def sep_distribution(theta_E, zs, vdf, cosmo=cosmo):
    theta_E_rad = theta_E / 206265
    paras = (zs,  vdf, cosmo)

    result = dblquad(sep_integral, 0, zs, 0, theta_E_rad, args=paras)

    return result[0]


'''
Functions for B
'''

def doublepowerlaw(M, phis, Ms, alpha, beta):

    phi =  phis / (10**(0.4*(M-Ms)*(alpha+1)) + 10**(0.4*(M-Ms)*(beta+1)))

    return phi

def Pmu_bright(mu):
    return 2/(mu-1)**3

def Pmu_total(mu):
    return 8 / mu**3

def N_Llim(Mlim, lumfun):
    '''
    Input Parameters:
        Llim: float
            The lower boundary of the luminosity
        lumfun: callable
            Call as Phi = lumfun(L)
    '''

    N = quad(lumfun, -40, Mlim)[0]
    return N

def magbias_differential(M, mu, lumfun, Pmu):
    return lumfun(M+2.5*np.log10(mu)) * Pmu(mu)


def magbias(Mlim, lumfun, Pmu, mu_min=2, mu_max = +np.inf):
    intfunc = lambda mu: N_Llim(Mlim+2.5*np.log10(mu), lumfun) * Pmu(mu)

    B = quad(intfunc, mu_min, mu_max)[0] / N_Llim(Mlim, lumfun)

    return B


'''
Functions for Number of Detectable Lenses
'''

def nqso_int(z, lumfun, Mlim):
    nqso = quad(lumfun, -40, Mlim)[0] * u.Mpc **-3

    dVdz = cosmo.differential_comoving_volume(z) * 4 * np.pi
    return nqso * dVdz


def Nqso(zmin, zmax, dz=0.01, mlim=22):
    totalnum = 0

    for z in np.arange(zmin, zmax+dz, dz):
        zmed = z + dz/2
        Phis = 10.9e-9 * 10**(-0.7*(zmed-6))
        alpha = -1.23
        beta = -2.73
        Ms = -24.9
        lumfun = lambda M: doublepowerlaw(M, Phis, Ms, alpha, beta)

        totalnum += nqso_int(zmed, lumfun, mlim) * dz

    return totalnum


def Nlens(zmin, zmax, dz=0.01, mlim=22):
    totalnum = 0

    for z in np.arange(zmin, zmax+dz, dz):
        zmed = z + dz/2
        Phis = 10.9e-9 * 10**(-0.7*(zmed-6))
        alpha = -1.23
        beta = -2.73
        Ms = -24.9
        lumfun = lambda M: doublepowerlaw(M, Phis, Ms, alpha, beta)

        Mlim = mlim - Mlim_VdB_accurate(zmed)
        tau = taum(zmed, vdf=vdf_default)
        B = magbias(Mlim, lumfun, Pmu=Pmu_bright)
        fm = Fmulti(tau, B)

        totalnum += nqso_int(zmed, lumfun, Mlim) * dz * fm

    return totalnum

def Fmulti(tau, B):
    return tau * B / (tau * B + 1 - tau)

def main():
    # examples of using these functions
    '''
    lensing optical depth
    '''
    # source redshift 
    ztest = 6

    # lensing optical depth
    tau = taum(ztest, vdf_analytical)
    print('Lensing Optical Depth at z=%.1f:'%ztest, tau)

    '''
    magnification bias
    '''
    # define the QLF
    Phis, Ms, alpha, beta = 10.9e-9, -25, -1.3, -2.6
    lumfun = lambda M: doublepowerlaw(M, Phis, Ms, alpha, beta)

    # Define the (absolute) magnitude limit
    absmaglim = -23

    # the magnification bias
    mb = magbias(absmaglim, lumfun, Pmu_bright)
    print('Magnification bias for Mabs=%.1f:'%absmaglim, mb)

    '''
    estimate number of detectable lensed quasars given a survey depth
    '''
    # the multiply-imaged fraction
    fm = Fmulti(tau, mb)

    # the number of detectable lensed quasar per redshift
    print('dNlens/dz (all sky) beyond Mabs=%.1f at z=%.1f'%(absmaglim, ztest),\
          nqso_int(ztest, lumfun, absmaglim).value * fm)

if __name__=='__main__':
    main()
