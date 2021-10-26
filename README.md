--- DOCUMENTATION FOR USING fd3_initiator -----
====================================================================================

fd3_initiator is a python wrapper around the code fd3 (http://sail.zpf.fer.hr/fdbinary/) for spectral disentangling.
The aim of this code was to provide a user-friendly way of using fd3. This code takes care 
of creating the master observation file and input files in the required format and runs fd3 
using these files.

Required packages and software:

1.GNU scientific Library for fd3

2.os [python]

3.glob [python]

4.numpy [python] 

To use fd3 in log of base 10:
gcc fd3fpolis_log10.c fd3sep.c triorb.c kepler.c mxfuns.c -lgsl -lm -lgslcblas -o fd3_log10
Use chmod +x fd3 (or fd3_log10) to make the code an executable and then use the main code (run_fd3.py).

--- HOW TO USE ----------
=========================

1. Copy the spectra from different epochs to the spectra/ folder named in order of epoch. We assume that your spectra file has two columns : wavelength (in log(Angstrom)) and normalised flux. 

2. The required parameters and settings for the disentangling process are initialised in paramfile.ini.

3. The RVs at each epoch and the light fraction (LF) are input as a table in the param_table.ini file.

4. Run "python run_fd3.py" in the terminal inside the folder.

5. Output files are generated named as the file identifier mentioned in paramfile.ini, with different 
extensions.

6
The disentangled spectra models (in rest frame) are in the .mod file. 
