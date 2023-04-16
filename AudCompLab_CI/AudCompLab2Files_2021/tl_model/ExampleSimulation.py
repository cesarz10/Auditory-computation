import numpy as np
from solve_one_cochlea import solve_one_cochlea
from cochlear_model2018 import *
import multiprocessing as mp
import hdf5storage

#%%
'''
CREATE THE STIMULI
'''
# create the stimulus
# Define variables
# create stimulus levels [0:10:100] as in matlab
L = np.arange(0., 101., 10.)
fs = 100e3
p0 = 2e-5
dur = 50e-3
t = np.arange(0., dur, 1./fs)
click_duration = 10 # 100 us click
stim = np.zeros((len(L), len(t)))
for j in range(len(L)):
    stim[j, 9 : 9 + click_duration] = 2 * np.sqrt(2) * p0 * 10**(L[j]/20)
   
#%%
'''
DEFINE ALL THE PARAMTERS HERE
'''
#Define all the variables for the model
sheraPdat = 'StartingPoles.dat'
poles = []
for line in open(sheraPdat, 'r'):
    poles.append(float(line.rstrip()))
sheraPo = np.array(poles)
probes = 'all' # 'half', 'all' or 'abr'
storeflag = 'evihmlbw'
probe_points=probes
Fs = fs
channels = np.min(stim.shape)
subjectNo = 1
sectionsNo = int(1e3)
output_folder = os.getcwd() + "/"
numH = 13.  
numM = 3.
numL = 3.
IrrPct = 0.05
nl = 'vel'
irregularities = 1
#%%
if __name__ == "__main__":
    # Create the opts file for passing the variables
    opts={}
    opts['sheraPo'] = sheraPo
    opts ['storeflag'] = storeflag
    opts ['probe_points'] = probe_points
    opts ['Fs'] = Fs
    opts ['channels'] = channels
    opts ['subjectNo'] = subjectNo
    opts ['sectionsNo'] = sectionsNo
    opts ['output_folder'] = output_folder
    opts ['numH'] = numH 
    opts ['numM'] = numM
    opts ['numL'] = numL
    opts ['IrrPct'] = IrrPct
    opts ['nl'] = nl
    opts ['L'] = L
    
    irr_on = irregularities * np.ones((1, channels)).astype('int')
    cochlear_list=[[cochlea_model(), stim[i], irr_on[0][i], i, opts] for i in range(channels)]
    
    print("running human auditory model 2018: Verhulst, Altoe, Vasilkov")
    p=mp.Pool(mp.cpu_count(),maxtasksperchild=1)
    output = p.map(solve_one_cochlea,cochlear_list)
    p.close()
    p.join()
    
    
    matcontent = {}
    matcontent[u'output'] = output
    destfilename = "Simulations.mat"
    destfile = output_folder + destfilename
    hdf5storage.write(matcontent, '.', destfile, store_python_metadata=True, matlab_compatible=True)

    print("cochlear simulation: done")
