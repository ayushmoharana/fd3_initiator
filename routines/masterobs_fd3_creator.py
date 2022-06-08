import numpy as np
import os, glob 


paramfile='paramfile.ini'
filefolder='spectra/'

os.chdir('../')
param=np.loadtxt(paramfile,unpack=True,dtype='str')
filename=str(param[2]+'.master.obs')
wv_start=float(param[0])
wv_end=float(param[1])

filelist=[]
wv=[]
flx=[]
print("")
print("---- Starting spectra collection: ", param[2] ,"  ----------")
os.chdir(filefolder)
#Look for spectra in the folder
for file in glob.glob("*.*"):
	filelist.append(file)
filelist=np.sort(filelist)

for i in range(len(filelist)):
	wex,fex=np.loadtxt(filelist[i],usecols=(0,1),unpack=True)
	wex_r=[]
	fex_r=[]
	for i in range(len(wex)):
		if (wex[i]<=wv_end) and (wex[i]>=wv_start):
			wex_r.append(wex[i])
			fex_r.append(fex[i])	
	wv.append(wex_r)
	flx.append(fex_r)
		
wv=np.asarray(wv)
flx=np.asarray(flx)	

input_array=[]

#Crop wavelength in the required range

for i in range(len(wv[0])):
	row=[]
	if (wv[0][i]<=wv_end) and (wv[0][i]>=wv_start):			 
		row.append(wv[0][i])
		for j in range(len(filelist)):
			try:
				row.append(float(flx[j][i]))
			except:
				print("--- ERROR : Wavelength range not available----")
				print("--- Cropping unavailable wavelegths ---")
				row=[]
	row=np.asarray(row)
			
#if wavelength present in all the spectra then select
	if len(row)==(len(filelist)+1):
		#print(row)
		input_array.append(row)

input_array=np.asarray(input_array)

#master file name
head=str(len(filelist)+1)+" X "+str(len(input_array))

os.chdir('../raw_fd3_output')
np.savetxt(filename,input_array,header=head,delimiter='	',fmt="%.8f")
print("")
print("-- Spectra list --")
for i in range(len(filelist)):	
	print(filelist[i])

print("")
print("---- Master observation file created -----")




