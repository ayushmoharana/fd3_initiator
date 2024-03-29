--- DOCUMENTATION FOR USING fd3_initiator -----
====================================================================================

fd3_initiator is a python wrapper around the code fd3 (http://sail.zpf.fer.hr/fdbinary/, cite: https://ui.adsabs.harvard.edu/abs/2004ASPC..318..111I/abstract) for spectral disentangling.
The aim of this code was to provide a user-friendly way of running fd3. This code takes care 
of creating the master observation file and input files in the required format and runs fd3 
using these files. Some corrections were added to the main fd3 executable to make it work for log of base 10 sampling.

This wrapper was used in:
1. Moharana et al., 2023, MNRAS (https://ui.adsabs.harvard.edu/abs/2023MNRAS.521.1908M/abstract)
2. Kahraman Aliçavuş et al.,2023, MNRAS (https://ui.adsabs.harvard.edu/abs/2023MNRAS.520.1601K/abstract)
3. Rozyczka et al., 2022, MNRAS (https://ui.adsabs.harvard.edu/abs/2022MNRAS.517.2485R/abstract)

Required packages and software:

1.GNU scientific Library for fd3

2.os [python]

3.glob [python]

4.numpy [python] 

To use fd3 in log of base 10, do this in modified_fd3_routines (skip this if the executable fd3_log10 works):
>gcc fd3fpolis_log10.c fd3sep.c triorb.c kepler.c mxfuns.c -lgsl -lm -lgslcblas -o fd3_log10

Copy fd3_log10 executable to raw_fd3_output.

Use chmod +x fd3 (or fd3_log10) to make the code an executable and then use the main code (run_fd3.py).

--- HOW TO USE ----------
=========================

1. Copy the spectra from different epochs to the spectra/ folder named in order of epochs. We assume that your spectra file has two columns : wavelength (in log(Angstrom) and uniformly sampled) and normalised flux. 

2. The required parameters and settings for the disentangling process are initialised in paramfile.ini.

3. The epoch of the spectra, heliocentric corrections, weight of the spectra and the light fraction (LF) are input as a table in the param_table.ini file.

4. Run "python run_fd3.py" in the terminal inside the folder.

5. Output files are generated named as the file identifier mentioned in paramfile.ini, with different 
extensions.

6. The disentangled spectra models (in rest frame) are in the .mod file and RVs in .rvs file. The output folder contains easy-to-use output files.
