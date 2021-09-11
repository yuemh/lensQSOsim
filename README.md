# lensQSOsim

This repository holds the code for generating the mock catalog of lensed quasars and the data products in Yue et al. (2021) (in prep).

## About the code

The code uses SIMQSO (https://github.com/imcgreer/simqso) to generate mock quasars, and depends on Numpy and Astropy. I use CosmoDC2 mock galaxy catalog (https://portal.nersc.gov/project/lsst/cosmoDC2/_README.html) to model the deflector population. I use GLAFIC (https://www.slac.stanford.edu/~oguri/glafic/, also see https://github.com/oguri/glafic2) to run lens modeling. I wrote a python interface for GLAFIC to run a large sample of lensing systems.

### If you would like to use the code, I strongly recommend contacting me before you try it. I am still improving this code to make it more flexible.

If you would like to run the code yourself, here is a brief instruction:

0. The CosmoDC2 mock galaxy catalog is too large to put in the repository. We thus include a python script get_cosmodc2_info.py to extract the required information from the whole CosmoDC2 catalog. You can download the whole CosmoDC2 catalog and run this script, or directly contact me at yuemh@email.arizona.edu.
1. Get all the dependencies, especially SIMQSO and GLAFIC. SIMQSO comes as a python package, while GLAFIC should be a binary executable file.
2. In setpaths.py, change the paths to the SIMQSO package and the executable GLAFIC file to where you download them.
4. Then you can start from mainrun.py, where the comments offer a good tutorial of generating a mock lensed quasar catalog with different settings.
5. The file qsopop.py contains the quasar luminosity functions (QLFs) we use. You can edit it to the QLFs you prefer.
6. The file lensing.py is the matching and lensing tools for the foreground galaxies and background quasars.
7. The file glafic.py is a python interface of the GLAFIC program. Don't edit it unless you know what you are doing! 


## About the data

In the lensQSOsim/data directory:

See the description file in lensQSOsim/data for a description of the columns.

The mock images are too large to put in this repository. Please contact me at yuemh@email.arizona.edu for them.
