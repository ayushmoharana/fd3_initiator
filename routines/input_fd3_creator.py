import numpy as np
import os, glob 

os.chdir('../')
paramfile='paramfile.ini'
rvlftablefile='param_table.ini'

param=np.loadtxt(paramfile,unpack=True,dtype='str')
try:
	t,rv,noise,lf_A,lf_B,lf_C=np.loadtxt(rvlftablefile,usecols=(0,1,2,3,4,5),dtype='str',unpack=True)
except:
	t,rv,noise,lf_A,lf_B=np.loadtxt(rvlftablefile,usecols=(0,1,2,3,4),dtype='str',unpack=True)

row_init=[]
matrix_rvlf=[]
row_trip=[]
row_bin=[]
row_opti=[]
row_init.append(str(param[2]+'.master.obs'))
inputfilename=str(param[2]+'_input.in')

param[0]=float(param[0])
param[1]=float(param[1])

for i in range(len(t)):
	try:		
		matrix_rvlf.append([t[i],rv[i],noise[i],lf_A[i],lf_B[i],lf_C[i]])
	except:
		matrix_rvlf.append([t[i],rv[i],noise[i],lf_A[i],lf_B[i]])

for i in range(len(param)):
	if i <=5 :
		row_init.append(param[i])

	if i >5 and i<=17: 
		row_trip.append(param[i])

	if i >17 and i<=31:
		row_bin.append(param[i])

	if i >31 and i<=34:
		row_opti.append(param[i])

os.chdir('raw_fd3_output')
row_init=np.asarray(row_init)
row_init=np.transpose(row_init)
f=open(inputfilename,'w')
np.savetxt(f, row_init,fmt='%s',newline='	')
f.write("\n")
for i in range(len(t)):
	f.write("\n")
	np.savetxt(f, matrix_rvlf[i], fmt='%s', newline='	')
f.write("\n")	
f.write("\n")
np.savetxt(f, row_trip, fmt='%s', newline=' ')
f.write("\n")
f.write("\n")
np.savetxt(f, row_bin, fmt='%s', newline=' ')
f.write("\n")
f.write("\n")
np.savetxt(f, row_opti, fmt='%s',newline=' ')
f.close()


print("")
print("---- fd3 input created -------------------")
