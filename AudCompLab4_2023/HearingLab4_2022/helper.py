# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 22:06:42 2020
@author: Tijmen Wartenberg
"""

import scipy.io.wavfile as wav
import numpy as np
import math
import os

def rms(insig, axis=None):
    """
    Calculates de root-mean-square of a waveform
    """
    return np.sqrt(np.mean(insig**2, axis=axis))

def Oct3smooth(vFin,vLin,vFout):
    """
    input:  vFin  [Hz]: frequency vector of the single-sided input spectrum
    vLin  [dB]: input spectrum
    vFout [Hz]: frequency vector of the smoothed spectrum
    """
    vLout = np.zeros(len(vFout))
    vLin = [10**(x/10) for x in vLin]
    
    for i in range(len(vLout)):                
      f1 = vFout[i]/(2**(1/6))             
      f2 = vFout[i]*(2**(1/6))             
      ind = np.where((vFin>=f1) & (vFin<=f2))
      indx = ind[0]
      vLout[i] = 10*math.log(np.mean(list( vLin[i] for i in indx )),10)
     
    return vLout

def SII(E, VE='normal',I=0,*args):
    """
    Calculate the Speech Intelligibility Index. Adapted to Python-code from the
    MATLAB script written by Hannes Muesch (2003).
        
    Parameters are passed to the procedure through pairs of "identifier" and corresponding "argument"
    Identifiesrs are always strings. Possible identifiers are:

     'E' Speech Spectrum Level (Section 3.6 in the standard)
     'VE' Vocal effort: {'normal','raised', 'loud', 'shout'} (by default: 'normal').
     'N' Equivalent Noise Spectrum Level (Section 3.15 in the standard)
     'T' Equivalent Hearing Threshold Level [dBHL] (Section 3.23 in the standard)
     'I' Band Importance function (Section 3.1 in the standard, integer between 0 and 4)
     'G' Insertion Gain [dB] (Section 3.28 in the standard)
    Except for 'E', which must be specified, all other parameters are optional.
    """
    try:
        N
    except NameError:
        N = -500+np.ones(18) # No noise is assumed
    except len(N) is not 18:
        print('Array should consist of eighteen numbers')
        
    try:
        T
    except NameError:
       T = np.zeros(18) #work
    except len(VE) is not 18:
        print('Array should consist of eighteen numbers')
        
    try:
        G
    except NameError:
       G = np.zeros(18) #work
    except len(G) is not 18:
        print('Array should consist of eighteen numbers')

    if not isinstance(VE,str):
        VE = 'normal' # default
        
    B = BandImportance(I)
    
    f = [160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000] #band center frequencies
    X = [0.6, -1.7, -3.9, -6.1, -8.2, -9.7, -10.8, -11.9, -12.5, -13.5, -15.4, -17.7, 
	-21.2, -24.2, -25.9, -23.6, -15.8, -7.1] #internal noise spectrum
    G = np.zeros(18) #work
    
    E = E + G
    V = E - 24
    B = np.maximum(V, N+G);
    C = 0.6*(B+([10*math.log10(x) for x in f]) -6.353)-80
    
    Z = np.zeros(18) 
    Z[0] = B[0]
    for i in range(1,18):
        y = [0.89*f[i]/x for x in f[0:(i)]]
        y = [math.log(x,10) for x in y]
        
        Z[i] = (10*math.log(10**(0.1*N[i]) + sum((10**(0.1*(B[0:(i)]+3.32*C[0:(i)]*y)))),10))
      
    
    L = 1 - (E - SpeechSpectrumLevel(VE) - 10)/160;
    L = (np.less_equal(L, np.ones(18))*1*L) + np.greater(L, np.ones(18))*1
    D = np.maximum(Z,X)

    K = (E-D+15)/30
    
    K = K*np.greater(K, np.zeros(18))*1
    K = (np.less_equal(K, np.ones(18))*1*K) + np.greater(K, np.ones(18))*1
   
    A = L*K
    S = np.sum(BandImportance(I)*A)
    print("The SII of the bands:", S)
    return S

def SpeechSpectrumLevel(VE):
    Ei = {'normal': [32.4100, 34.4800, 34.7500, 33.9800, 34.5900, 34.2700, 32.0600, 28.3000, 25.0100, 23.0000, 20.1500, 17.3200, 13.1800, 11.5500, 9.3300, 5.3100, 2.5900, 1.1300],
          'raised': [33.8100, 33.9200, 38.9800, 38.5700, 39.1100, 40.1500, 38.7800, 36.3700, 33.8600, 31.8900, 28.5800, 25.3200, 22.3500, 20.1500, 16.7800, 11.4700, 7.6700, 5.0700],
          'loud':  [35.2900, 37.7600, 41.5500, 43.7800, 43.3000, 44.8500, 45.5500, 44.0500, 42.1600, 40.5300, 37.7000, 34.3900, 30.9800, 28.2100, 25.4100, 18.3500, 13.8700, 11.3900],
         'shout':  [30.7700, 36.6500, 42.5000, 46.5100, 47.4000, 49.2400, 51.2100, 51.4400, 51.3100, 49.6300, 47.6500, 44.3200, 40.8000, 38.1300, 34.4100, 28.2400, 23.4500, 20.7200]}
    E = Ei[VE]
    return E
    
def BandImportance(I):
        import numpy
        BIarr = numpy.array([[0.0083, 0.    , 0.0365, 0.0168, 0.    , 0.0114, 0.    ],
        [0.0095, 0.    , 0.0279, 0.013 , 0.024 , 0.0153, 0.0255],
        [0.015 , 0.0153, 0.0405, 0.0211, 0.033 , 0.0179, 0.0256],
        [0.0289, 0.0284, 0.05  , 0.0344, 0.039 , 0.0558, 0.036 ],
        [0.044 , 0.0363, 0.053 , 0.0517, 0.0571, 0.0898, 0.0362],
        [0.0578, 0.0422, 0.0518, 0.0737, 0.0691, 0.0944, 0.0514],
        [0.0653, 0.0509, 0.0514, 0.0658, 0.0781, 0.0709, 0.0616],
        [0.0711, 0.0584, 0.0575, 0.0644, 0.0751, 0.066 , 0.077 ],
        [0.0818, 0.0667, 0.0717, 0.0664, 0.0781, 0.0628, 0.0718],
        [0.0844, 0.0774, 0.0873, 0.0802, 0.0811, 0.0672, 0.0718],
        [0.0882, 0.0893, 0.0902, 0.0987, 0.0961, 0.0747, 0.1075],
        [0.0898, 0.1104, 0.0938, 0.1171, 0.0901, 0.0755, 0.0921],
        [0.0868, 0.112 , 0.0928, 0.0932, 0.0781, 0.082 , 0.1026],
        [0.0844, 0.0981, 0.0678, 0.0783, 0.0691, 0.0808, 0.0922],
        [0.0771, 0.0867, 0.0498, 0.0562, 0.048 , 0.0483, 0.0719],
        [0.0527, 0.0728, 0.0312, 0.0337, 0.033 , 0.0453, 0.0461],
        [0.0364, 0.0551, 0.0215, 0.0177, 0.027 , 0.0274, 0.0306],
        [0.0185, 0.    , 0.0253, 0.0176, 0.024 , 0.0145, 0.0    ]])
        B = [BIarr[:,I]]
        return B

def audioread(fname):
    """
    Reads an audio file, scaling the output data to values between -1 and 1
    """
    # import scipy.io.wavfile as wav

    fs,insig = wav.read(fname)
    if insig.dtype == 'int16':
        nb_bits = 16 # -> 16-bit wav files
    elif insig.dtype == 'int32':
        nb_bits = 32 # -> 32-bit wav files
    max_nb_bit = float(2 ** (nb_bits - 1))
    insig = insig / (max_nb_bit + 1.0)
    return fs, insig

def SpeechArray(file_path,nr_of_sentences=10):
    """
    file_path = location of wav files (full path)
    nr_of_sentences = number of sentences to be read, by default 10
    """
    print(file_path)
    lijsten = os.listdir(file_path)
    wavarray = []
    
    print('\t {:.0f} sounds were found'.format(len(lijsten)))
    for l in lijsten:
        if 'noise' in l: # Skip noise file
            print('Noise file skipped {}'.format(l))
            continue
        else:
            print('Loading sentence {}'.format(l))
        fs,sig = audioread(file_path+l)
        wavarray.append([sig])
    
    print('\t A subset of {:.0f} sounds will be returned'.format(nr_of_sentences))    
    sentences = np.concatenate(np.random.randint(len(wavarray), size=(1,nr_of_sentences)), axis=0) # choose twenty random lists

    outsigs = (np.concatenate([wavarray[s] for s in sentences],axis=1))
    outsigs = outsigs.transpose()
    outsigs = outsigs[:,0]
    return outsigs

def audiowrite(insig,fs,fname):
    # wav.write(fname, fs, np.array(insig,dtype = np.int16))
    nb_bits = 16 # -> 16-bit wav files
    max_nb_bit = float(2 ** (nb_bits - 1))
            
    insig2store = np.array(insig*max_nb_bit,dtype = np.int16)
    wav.write(fname, fs, insig2store)