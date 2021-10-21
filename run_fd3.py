import os
import numpy as np

paramfile='paramfile.ini'
param=np.loadtxt(paramfile,unpack=True,dtype='str')
#for i in range(len(param)):	
#	print(i,param[i])
os.chdir('routines/')
inputfilename=str(param[2]+'_input.in')
outputfilename=str(param[2]+'_output.out')

run_masterfile=str("python "+"masterobs_fd3_creator.py")
os.system(run_masterfile)

run_inputgen=str("python "+"input_fd3_creator.py")
os.system(run_inputgen)

os.chdir('../raw_fd3_output')
if param[35]=='e':
	if param[36]=='Y':	
		runfd3=str("./fd3 < "+inputfilename +" > "+outputfilename)
	if param[36]=='N':
		runfd3=str("./fd3 < "+inputfilename)#+" > "+outputfilename)
else:
	if param[36]=='Y':
		runfd3=str("./fd3_log10 < "+inputfilename +" > "+outputfilename)
	if param[36]=='N':
		runfd3=str("./fd3_log10 < "+inputfilename)#+" > "+outputfilename)
os.system(runfd3)

os.chdir('../routines/')
create_out=str("python "+"output_creator.py")
os.system(create_out)
