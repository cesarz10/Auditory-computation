# python clone of ExampleAnalysis.m
# Deepak Baby, UGent, Aug 2018

import numpy as np
import hdf5storage
#import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

def rms (x):
    # compute rms of a matrix
    sq = np.mean(np.square(x), axis = 0)
    return np.sqrt(sq)

##############################################
# stimulus related stuff
L = np.arange(0., 101., 10.)
p0 = 2e-5
No = 245 # pick a channel number
# load Simulation results
# we will store it in a dictionary
simulations = hdf5storage.loadmat("Simulations.mat")
output = simulations['output']
#%%    
CF = output[0]['cf'].T # we will take transpose of it for plotting
plt.figure(1)
plt.plot(CF)
plt.xlabel('Cochlear Channel Number [-]')
plt.ylabel('Characteristic Frequency [Hz]')
plt.title('Characteristic frequencies')
plt.grid()

fs_c = output[0]['fs_bm']
fs=output[0]['fs_an']
t_c = np.arange(output[0]['v'].shape[0]) / fs_c
t = np.arange(output[0]['anfH'].shape[0]) / fs
f_c = np.arange(0, fs_c, fs_c/output[0]['v'].shape[0])
f = np.arange(0, fs, fs/output[0]['anfH'].shape[0])

# reorganise data for easier processing
num_cf = CF.shape[0]
numsamples = len(t_c)
num_levels = len(L)

oae = np.zeros((numsamples, num_levels))
vrms = np.zeros((num_cf, num_levels))
ihcrms = np.zeros((num_cf, num_levels))
v = np.zeros((numsamples, num_levels))
ihc = np.zeros((numsamples, num_levels))

for i in range(len(L)):
    oae[:,i] = output[i]['e']
    vrms[:,i] = rms(output[i]['v'])
    ihcrms [:,i] = rms(output[i]['ihc'])
    v[:,i] = output[i]['v'][:,No]
    ihc[:,i] = output[i]['ihc'][:,No]
    
#%% the OAE
#this figure is the ear-canal pressure: 
#to get the reflection-source OAE, do the following simulation:
#OAE_{reflections on}-OAE_{reflections off}
#to get the distortion-source OAE, do a normalisation using a low level
#(linear) simulation that the reflections off
legend_labels=[]
for i in range(num_levels) :
    legend_labels.append(str(int(L[i])))
    
plt.figure(2)
plt.subplot(2,1,1)
lineObjects = plt.plot(1000*t_c, oae[:,:])
plt.xlabel('Time [ms]')
plt.ylabel('Ear Canal Pressure [Pa]')
plt.xlim((0, 20)), plt.ylim((-0.02, 0.02))
plt.title('OAE')
plt.legend(lineObjects, legend_labels)
plt.grid()

# zero pad oae to plot spectrum
oae[:200,:] = 0
plt.subplot(2,1,2)
lineObjects = plt.plot(f_c/1000, 20*np.log10(np.abs(np.fft.fft(oae.T/p0)).T));
plt.xlabel('Frequency [kHz]')
plt.ylabel('EC Magnitude [dB re p0]')
plt.legend(lineObjects, legend_labels)
plt.xlim((0, 12))
plt.grid()

#%% v_bm and V_IHC
plt.figure(3)
ax1 = plt.subplot(2, 2, 1)
lineObj = ax1.plot(1000*t_c,v)
plt.xlabel('Time [ms]'), plt.ylabel('v_{bm} [m/s]'), plt.xlim((0,30))
plt.legend(lineObj, legend_labels)
ax1.set_title("CF of " + str(CF[No])+ " Hz")
ax2 = plt.subplot(2,2,2)
lineObj = ax2.plot(CF/1000,20*np.log10(vrms)),
plt.xlabel('CF [kHz]'), plt.ylabel('rms of v_{bm} [dB re 1 m/s]'),
plt.xlim((0, 14))
plt.ylim((np.max(np.max(20*np.log10(vrms)))-100, np.max(np.max(20*np.log10(vrms)))+10))
ax2.set_title('Excitation Pattern')
ax3 = plt.subplot(2,2,3)
lineObj = ax3.plot(1000*t_c, ihc)
plt.xlabel('Time [ms]'), plt.ylabel('V_{ihc} [V]')
plt.xlim((0, 30))
plt.legend(lineObj, legend_labels)
ax3.set_title("CF of '" + str(CF[No]) + " Hz")
ax4 = plt.subplot(2,2,4)
ax4.plot(CF/1000, ihcrms)
plt.xlabel('CF [kHz]'), plt.ylabel('rms of V_{ihc} [V]'),
plt.xlim((0,14))
#plt.ylim((np.max(np.max(20*np.log10(ihcrms)))-100, np.max(np.max(20*np.log10(ihcrms)))+10))
ax4.set_title('Excitation Pattern')
#%%
# Reorganising for easier processing
num_samples_resamp =  output[0]['anfH'].shape[0] # number of downsampled signal points
HSR = np.zeros((num_samples_resamp, num_levels))
MSR = np.zeros((num_samples_resamp, num_levels))
LSR = np.zeros((num_samples_resamp, num_levels))
AN = np.zeros((num_samples_resamp, num_levels))
CN = np.zeros((num_samples_resamp, num_levels))
IC = np.zeros((num_samples_resamp, num_levels))
W1 = np.zeros((num_samples_resamp, num_levels))
W3 = np.zeros((num_samples_resamp, num_levels))
W5 = np.zeros((num_samples_resamp, num_levels))
EFR = np.zeros((num_samples_resamp, num_levels))

for n in range(num_levels):
    HSR[:,n] = output[n]['anfH'][:,No]
    MSR[:,n] = output[n]['anfM'][:,No]
    LSR[:,n] = output[n]['anfL'][:,No]
    AN[:,n]= output[n]['an_summed'][:,No]
    CN[:,n] = output[n]['cn'][:,No]
    IC[:,n] = output[n]['ic'][:,No]
    W1[:,n] = output[n]['w1']
    W3[:,n] = output[n]['w3']
    W5[:,n] = output[n]['w5']
    EFR[:,n] = output[n]['w1'] + output[n]['w3'] + output[n]['w5']


# single unit responses
plt.figure(4)
ax1 = plt.subplot(3,2,1)
lineObj = ax1.plot(1000*t, HSR)
ax1.set_title("CF of " + str(CF[No]) + " Hz")
plt.xlim((0, 20)), plt.xlabel("Time [ms]"), plt.ylabel('HSR fiber [spikes/s]')
plt.legend(lineObj, legend_labels)
plt.subplot(3,2,3), plt.plot(1000*t, MSR)
plt.xlim((0, 20)), plt.xlabel('Time [ms]'), plt.ylabel('HSR fiber [spikes/s]')
plt.subplot(3,2,5), plt.plot(1000*t, LSR)
plt.xlim((0, 20)), plt.xlabel('Time [ms]'), plt.ylabel('HSR fiber [spikes/s]')

ax2 = plt.subplot(3,2,2)
ax2.plot(1000*t, AN)
ax2.set_title("CF of " + str(CF[No]) + " Hz")
plt.xlim((0, 20)), plt.xlabel('Time [ms]'), plt.ylabel('sum AN [spikes/s]')
#Spikes summed across all fibers @ 1 CF
plt.subplot(3,2,4), plt.plot(1000*t, CN)
plt.xlim((0, 20)),plt.xlabel('Time [ms]'),plt.ylabel('CN [spikes/s]')
plt.subplot(3,2,6), plt.plot(1000*t, IC)
plt.xlim((0,20)), plt.xlabel('Time [ms]'),plt.ylabel('IC [spikes/s]')


#population responses
plt.figure(5)
ax1 = plt.subplot(4,1,1)
lineObj = ax1.plot(1000*t, 1e6*W1)
plt.xlim((0, 20))
plt.xlabel('Time [ms]'), plt.ylabel('W-1 [\muV]')
plt.legend(lineObj, legend_labels)
ax1.set_title("Population Responses summed across simulated CFs")
plt.subplot(4,1,2), plt.plot(1000*t, 1e6*W3)
plt.xlim((0,20)), plt.xlabel('Time [ms]'), plt.ylabel('W-3 [\muV]')
plt.subplot(4,1,3), plt.plot(1000*t, 1e6*W5)
plt.xlim((0, 20)), plt.xlabel('Time [ms]'), plt.ylabel('W-5 [\muV]')
plt.subplot(4,1,4), plt.plot(1000*t, 1e6*EFR)
plt.xlim((0, 20)), plt.xlabel('Time [ms]'), plt.ylabel('EFR [\muV]')

plt.show()


#The model code and interface was written by Alessandro Altoe and Sarah Verhulst (copyright 2012,2014,2015,2016,2018) 
#and is licensed under the UGent acadamic license (see details in license file that is part of this repository). 
#The Verhulstetal2018Model consists of the following files: 
#tridiag.so, cochlea_utils.c, run_model2018.py, model2018.m, cochlear_model2017.py, inner_hair_cell2018.py, auditory_nerve2017.py, ic_cn2017.py, ExampleSimulation.m, ExampleAnalysis.m, the HI profiles in the Poles folder. 
 
# This is a python clone written by Deepak Baby for the interface code written in matlab
