# The mock catalog of lensed quasars at z_qso<7.5

lowzlensedqso.fits: the mock catalog for quasars at z_qso<5. Corresponding to one all-sky area.

highzlensedqso.fits: the mock catalog for quasars at z_qso>5. Corresponding to 20 all-sky area.

qsospec.fits: the simulated (un-magnified) spectrum of the mock quasars.

## Columns in the catalog

### Quasar-related 

qra, qdec: the coordinate of the background quasar.

zq: the redshift of the quasar.

absMag: the absolute magnitude at rest-frame 1450 \AA of the quasar (not magnified).

qidx: unique identifier of the quasar, which is used to get its spectrum.

### Galaxy-related (see CosmoDC2 documents for more information)

gra, gdec: the coordinate of the foreground deflector galaxy, in deg.

zg: the redshift of the galaxy

size_true, size_minor_true: the semi-major and semi-minor axis of the galaxy, in arcsec

pa_lens: the position angle of the galaxy (follows GLAFIC conversion)

stellar_mass: the stellar mass in M_sun of the galaxy

sigma: the velocity dispersion

galaxy_id: the unique identifier of the galaxy in the CosmoDC2 catalog


### Lensing-related

convergence: the external convergence

shear_1, shear_2: the x- and y- component of the external shear

