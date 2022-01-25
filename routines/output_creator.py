import numpy as np
import os, glob 
import math


paramfile='paramfile.ini'
rvlftablefile='param_table.ini'
filefolder='raw_fd3_output/'

os.chdir('../')
param=np.loadtxt(paramfile,unpack=True,dtype='str')
t=np.loadtxt(rvlftablefile,usecols=(0),unpack=True)

log_base=param[35]
modelfile=str(param[2]+'.mod')
rvfile=str(param[2]+'.rvs')

os.chdir(filefolder)
#Output Model file
if int(param[5])==0:
	logwv,flx1,flx2=np.loadtxt(modelfile,usecols=(0,1,2),unpack=True)
if int(param[5])==1:
	logwv,flx1,flx2,flx3=np.loadtxt(modelfile,usecols=(0,1,2,3),unpack=True)

#Output RV file
if int(param[5])==0:
	rv1,rv2=np.loadtxt(rvfile,usecols=(0,1),unpack=True)
if int(param[5])==1:
	rv1,rv2,rv3=np.loadtxt(rvfile,usecols=(0,1,2),unpack=True)
	
#save new output
mod_speclist=[]
rvlist=[]

for i in range(len(logwv)):
	if log_base=='e':
		if int(param[5])==0:
			mod_speclist.append([math.exp(logwv[i]),flx1[i],flx2[i]])
		if int(param[5])==1:
			mod_speclist.append([math.exp(logwv[i]),flx1[i],flx2[i],flx3[i]])
	else:
		if int(param[5])==0:
			mod_speclist.append([pow(10,logwv[i]),flx1[i],flx2[i]])
		if int(param[5])==1:
			mod_speclist.append([pow(10,logwv[i]),flx1[i],flx2[i],flx3[i]])


for i in range(len(t)):
	if int(param[5])==0:
		rvlist.append([t[i],rv1[i],rv2[i]])
	if int(param[5])==1:
		rvlist.append([t[i],rv1[i],rv2[i],rv3[i]])
		
os.chdir('../output/')
if int(param[5])==0:
	np.savetxt(modelfile+'.out',mod_speclist,header="WV_ang	Flux_A	Flux_B",delimiter='	',fmt="%.8f")
if int(param[5])==1:
	np.savetxt(modelfile+'.out',mod_speclist,header="WV_ang	Flux_A	Flux_B	Flux_C",delimiter='	',fmt="%.8f")

if int(param[5])==0:
	np.savetxt(rvfile+'.out',rvlist,header="time	RV_A	RV_B",delimiter='	',fmt="%.8f")
if int(param[5])==1:
	np.savetxt(rvfile+'.out',rvlist,header="time	RV_A	RV_B	RV_C",delimiter='	',fmt="%.8f")

with open(rvfile+'.out',"a+") as file_object:
	file_object.write("\n")
	file_object.write("# Multiply the RV step per bin (check output on the terminal or the file) to the RVs to get physical values")
