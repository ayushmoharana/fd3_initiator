import numpy as np
import os, glob 
import math


paramfile='paramfile.ini'
rvlftablefile='param_table.ini'
filefolder='raw_fd3_output/'

def error_from_res2(res_table,lfA,lfB):
    error_A=[]
    error_B=[]
    error_C=[]
    wv=[]
    for i in range(len(res_table)):
        sumA=0
        sumB=0
        sumC=0
        for j in range(len(lfA)):
            sumA=sumA+(abs(np.asarray(res_table[i][j+1]))/np.sqrt(lfA[j]))
            sumB=sumB+(abs(np.asarray(res_table[i][j+1]))/np.sqrt(lfB[j]))
        error_A.append(sumA/len(lfA))
        error_B.append(sumB/len(lfA))
        #wv.append(res_table[i][0])

    return error_A,error_B

def error_from_res3(res_table,lfA,lfB,lfC):
    error_A=[]
    error_B=[]
    error_C=[]
    wv=[]
    for i in range(len(res_table)):
        sumA=0
        sumB=0
        sumC=0
        for j in range(len(lfA)):
            sumA=sumA+(abs(np.asarray(res_table[i][j+1]))/np.sqrt(lfA[j]))
            sumB=sumB+(abs(np.asarray(res_table[i][j+1]))/np.sqrt(lfB[j]))
            sumC=sumC+(abs(np.asarray(res_table[i][j+1]))/np.sqrt(lfC[j]))
        error_A.append(sumA/len(lfA))
        error_B.append(sumB/len(lfA))
        error_C.append(sumC/len(lfA))
        #wv.append(res_table[i][0])

    return error_A,error_B,error_C

os.chdir('../')
param=np.loadtxt(paramfile,unpack=True,dtype='str')
t=np.loadtxt(rvlftablefile,usecols=(0),unpack=True)

if int(param[5])==0:
	lfA,lfB=np.loadtxt(rvlftablefile,usecols=(3,4),unpack=True)
if int(param[5])==1:
	lfA,lfB,lfC=np.loadtxt(rvlftablefile,usecols=(3,4,5),unpack=True)

log_base=param[35]
modelfile=str(param[2]+'.mod')
rvfile=str(param[2]+'.rvs')
residual_file=str(param[2]+'.res')

os.chdir(filefolder)
#Output Model file
if int(param[5])==0:
	logwv,flx1,flx2=np.loadtxt(modelfile,usecols=(0,1,2),unpack=True)
	res_table=np.loadtxt(residual_file)
	er1,er2=error_from_res2(res_table,lfA,lfB)
if int(param[5])==1:
	logwv,flx1,flx2,flx3=np.loadtxt(modelfile,usecols=(0,1,2,3),unpack=True)
	res_table=np.loadtxt(residual_file)
	er1,er2,er3=error_from_res3(res_table,lfA,lfB,lfC)

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
			mod_speclist.append([math.exp(logwv[i]),flx1[i],er1[i],flx2[i],er2[i]])
		if int(param[5])==1:
			mod_speclist.append([math.exp(logwv[i]),flx1[i],er1[i],flx2[i],er2[i],flx3[i],er3[i]])
	else:
		if int(param[5])==0:
			mod_speclist.append([pow(10,logwv[i]),flx1[i],er1[i],flx2[i],er2[i]])
		if int(param[5])==1:
			mod_speclist.append([pow(10,logwv[i]),flx1[i],er1[i],flx2[i],er2[i],flx3[i],er3[i]])


for i in range(len(t)):
	if int(param[5])==0:
		rvlist.append([t[i],rv1[i],rv2[i]])
	if int(param[5])==1:
		rvlist.append([t[i],rv1[i],rv2[i],rv3[i]])
		
os.chdir('../output/')
if int(param[5])==0:
	np.savetxt(modelfile+'.out',mod_speclist,header="WV_ang    Flux_A    Fl_A_er    Flux_B    Fl_B_er",delimiter='	',fmt="%.8f")
if int(param[5])==1:
	np.savetxt(modelfile+'.out',mod_speclist,header="WV_ang    Flux_A    Fl_A_er    Flux_B    Fl_B_er   Flux_C   Fl_C_er",delimiter='	',fmt="%.8f")

if int(param[5])==0:
	np.savetxt(rvfile+'.out',rvlist,header="time	RV_A	RV_B",delimiter='	',fmt="%.8f")
if int(param[5])==1:
	np.savetxt(rvfile+'.out',rvlist,header="time	RV_A	RV_B	RV_C",delimiter='	',fmt="%.8f")

with open(rvfile+'.out',"a+") as file_object:
	file_object.write("\n")
	file_object.write("# Multiply the RV step per bin (check output on the terminal or the file) to the RVs to get physical values")
