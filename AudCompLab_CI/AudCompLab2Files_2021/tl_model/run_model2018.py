import numpy as np
import scipy as sp
import scipy.io as sio
from cochlear_model2018 import *
import os
import warnings
import multiprocessing as mp
import ctypes as c
import time
import sys
import inner_hair_cell2018 as ihc
import auditory_nerve2018 as anf
import ic_cn2018 as nuclei

#this relates to python 3.6 on ubuntu
#there is one future warning related to "scipy.signal.decimate" in this file
#there is one runtime warning related to firwin "scipy.signal.decimate" in ic_cn2017.py (not important)
#so we suppress these warnings here
warnings.filterwarnings("ignore")

par=sio.loadmat('input.mat')
probes=np.array(par['probes']) 
storeflag_in=np.array(par['storeflag'],dtype=str)
storeflag=storeflag_in[0]
probe_points=probes
Fs=par['Fs']
Fs=Fs[0][0]
stim=par['stim']
channels=par['channels']
channels=channels[0][0]
subjectNo=int(par['subject'])
sectionsNo=int(par['sectionsNo'])
t_f=(par['data_folder'])
output_folder=str(t_f[0])
lgt=len(stim[0])
sheraPo=par['sheraPo']
if(max(np.shape(sheraPo))==1):
    sheraPo=sheraPo[0][0]
else:
    sheraPo=sheraPo[:,0]
numH=par['nH']
if(max(np.shape(numH))==1):
    numH=numH[0][0]
else:
    numH=numH[:,0]
numM=par['nM']
if(max(np.shape(numM))==1):
    numM=numM[0][0]
else:
    numM=numM[:,0]
numL=par['nL']
if(max(np.shape(numL))==1):
    numL=numL[0][0]
else:
    numL=numL[:,0]




IrrPct=par['IrrPct']
IrrPct=IrrPct[0][0]
nl=np.array(par['non_linear_type'])
#print(IrrPct)
#print(sheraPo)
irr_on=np.array(par['irregularities'])
d=len(stim[0].transpose())
print("running human auditory model 2018: Verhulst, Altoe, Vasilkov")
sig=stim

cochlear_list=[ [cochlea_model(),sig[i],irr_on[0][i],i] for i in range(channels)]
#sheraPo = np.loadtxt('StartingPoles.dat', delimiter=',')
#print(sheraPo)
#%%
def solve_one_cochlea(model): #definition here, to have all the parameter implicit
    ii=model[3]
    coch=model[0]
    sig=model[1]
    coch.init_model(model[1],Fs,sectionsNo,probe_points,Zweig_irregularities=model[2],sheraPo=sheraPo,subject=subjectNo,IrrPct=IrrPct,non_linearity_type=nl)

    coch.solve()
    magic_constant=0.118;
    Vm=ihc.inner_hair_cell_potential(coch.Vsolution*magic_constant,Fs)
    dec_factor=5
    Vm_resampled=sp.signal.decimate(Vm,dec_factor,axis=0,n=30,ftype='fir')
    Vm_resampled[0:5,:]=Vm[0,0]; #resting value to eliminate noise from decimate
    Fs_res=Fs/dec_factor
#    print(np.shape(coch.Vsolution),np.shape(Vm_resampled))
    if 'v' in storeflag:
        f=open(output_folder+"v"+str(ii+1)+".np",'w')
        coch.Vsolution.tofile(f)
        f.close()
    if 'y' in storeflag:
        f=open(output_folder+"y"+str(ii+1)+".np",'w')
        coch.Ysolution.tofile(f)
        f.close()
    if 'i' in storeflag:
        f=open(output_folder+"ihc"+str(ii+1)+".np",'w')
        Vm.tofile(f)
        f.close()
    if 'h' in storeflag or 'b' in storeflag:
        anfH=anf.auditory_nerve_fiber(Vm_resampled,Fs_res,2)*Fs_res
#print(np.shape(coch.Vsolution),np.shape(Vm_resampled),np.shape(anfH))

    if 'h' in storeflag:
        f=open(output_folder+"anfH"+str(ii+1)+".np",'w')
        anfH.tofile(f)
        f.close()
    if 'm' in storeflag or 'b' in storeflag:
        anfM=anf.auditory_nerve_fiber(Vm_resampled,Fs_res,1)*Fs_res
    if 'm' in storeflag:
        f=open(output_folder+"anfM"+str(ii+1)+".np",'w')
        anfM.tofile(f)
        f.close()
    if 'l' in storeflag or 'b' in storeflag:
        anfL=anf.auditory_nerve_fiber(Vm_resampled,Fs_res,0)*Fs_res
    if 'l' in storeflag:
        f=open(output_folder+"anfL"+str(ii+1)+".np",'w')
        anfL.tofile(f)
        f.close()
    if 'e' in storeflag:
        f=open(output_folder+"emission"+str(ii+1)+".np",'w')
        coch.oto_emission.tofile(f)
        f.close()
    if 'b' in storeflag or 'w' in storeflag:
        cn,anSummed=nuclei.cochlearNuclei(anfH,anfM,anfL,numH,numM,numL,Fs_res)
        ic=nuclei.inferiorColliculus(cn,Fs_res)
        if 'b' in storeflag:
            f=open(output_folder+"cn"+str(ii+1)+".np",'w')
            cn.tofile(f)
            f.close()
            f=open(output_folder+"AN"+str(ii+1)+".np",'w')
            anSummed.tofile(f)
            f.close()
            f=open(output_folder+"ic"+str(ii+1)+".np",'w')
            ic.tofile(f)
            f.close()
        if 'w' in storeflag:
            w1=nuclei.M1*np.sum(anSummed,axis=1);
            w3=nuclei.M3*np.sum(cn,axis=1)
            w5=nuclei.M5*np.sum(ic,axis=1)
            f=open(output_folder+"1w"+str(ii+1)+".np",'w')
            w1.tofile(f)
            f.close()
            f=open(output_folder+"3w"+str(ii+1)+".np",'w')
            w3.tofile(f)
            f.close()
            f=open(output_folder+"5w"+str(ii+1)+".np",'w')
            w5.tofile(f)
            f.close()
    f=open(output_folder+"cf"+str(ii+1)+".np",'w')
    coch.cf.tofile(f)
    f.close()



if __name__ == "__main__":
    s1=time.clock()
    p=mp.Pool(mp.cpu_count(),maxtasksperchild=1)
    p.map(solve_one_cochlea,cochlear_list)
    p.close()
    p.join()

    print("cochlear simulation: done")
