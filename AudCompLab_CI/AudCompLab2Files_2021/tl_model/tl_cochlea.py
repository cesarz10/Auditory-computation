"""
Generate transmission line cochlear outputs
"""

import numpy as np
from .cochlear_model2018 import *

import os

def tl_cochlea(x):
    '''
    Generates the transmission cochlear model outputs for stimuli x
    '''
    # DEFINE ALL THE PARAMTERS HERE
    
    #Define all the variables for the model
    basedir = os.getcwd()
    #sheraPdat = os.path.join(basedir, 'tl_model/Poles/Flat35/StartingPoles.dat')
    #sheraPdat = os.path.join(basedir, 'tl_model/Poles/Slope35/StartingPoles.dat')
    sheraPdat = os.path.join(basedir, 'tl_model/StartingPoles.dat')
    poles = []
    for line in open(sheraPdat, 'r'):
        poles.append(float(line.rstrip()))
    sheraPo = np.array(poles)
    # Create the opts file for passing the variables
    probe_points = 'abr'
    Fs = 100e3
    sectionsNo = int(1e3)
    IrrPct = 0.05
    nl = 'vel'
    
    print ("Running human auditory model 2018: Verhulst, Altoe, Vasilkov")
    coch = cochlea_model()
    coch.init_model(x, Fs,sectionsNo, probe_points, Zweig_irregularities=1,sheraPo=sheraPo, subject=1, IrrPct=IrrPct, non_linearity_type=nl)
    coch.solve()
    print("cochlear simulation: done")
    
    output = {}
    output['cf'] = coch.cf
    output['fs'] = Fs
    output['bmm'] = coch.Vsolution
    
    return output
