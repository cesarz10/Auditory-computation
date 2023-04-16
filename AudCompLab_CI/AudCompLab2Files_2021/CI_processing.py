#!/usr/bin/env python

"""
Paradigm of the processing strategy applied in Cochlear Implants. 
A monophasic or a biphasic pulse is used to excite the CI electrodes 
of the cochlear implant. The pulse is modulated according to the 
acoustic amplitude of the input signal.
Adapted from the High Resolution Strategy (HiRes).

Fotis Drakopoulos, UGent, March 2019

signal, fsy: input signal to be processed and its sampling frequency
(default = 16 kHz)- the signal will be downsampled to 16 kHz for the
processing.

M: the size of the filterbank (or the amount of the respective 
electrodes) - if 1 is given then no fb is used
Ts: duration of a stimulation cycle in secs - default is 1 ms

MCL, THL, IDR: If these arguments are provided then a logarithmic 
compression function is used which ensures that the envelope outputs fit 
the patient's dynamic range. The compression is applied according to the 
HiRes strategy. If no arguments are given then no compression is applied.
MCL, THL: Most Comfortable Level and Threshold in dB-SPL for each 
frequency band - the arrays must have the same size as the value of M. 
IDR: The input dynamic range (set by the clinicial) - default is 60 dB.

Modified by Fotis Drakopoulos, UGent, March 2020
"""

from __future__ import division, print_function
import sys
from argparse import ArgumentParser
import numpy as np
import scipy.io.wavfile
from scipy.signal import square, butter, filtfilt, resample_poly

def build_argparser():
    parser = ArgumentParser()
    parser.add_argument("-i", "--wavfile", help="Input wavefile to be processed", required=True, type=str)
    parser.add_argument("-m", "--bands", help="Size of the filterbank", default=16, type=int)
    parser.add_argument("-ts", "--pulse_duration", help="Duration of a stimulation cycle in secs", default=0.001, type=float)
    parser.add_argument("-p", "--pulse", help="1 for monophasic pulse, 2 for biphasic pulse", default=2, type=int)
    parser.add_argument("-mcl", "--comfortable_levels", nargs='+', help="List containing the Most Comfortable Level in dB-SPL for each frequency band", default=[])
    parser.add_argument("-thl", "--threshold", nargs='+', help="List containing the Hearing Thresholds in dB-SPL for each frequency band", default=[])
    parser.add_argument("-idr", "--dynamic_range", help="Input Dynamic Range", default=60, type=float)

    return parser

def wavfile_read(wavfile,fs=[]):
    # if fs is given the signal is resampled to the given sampling frequency
    fs_signal, signal = scipy.io.wavfile.read(wavfile)
    if not fs:
        fs=fs_signal

    if signal.dtype == 'int16':
        nb_bits = 16 # -> 16-bit wav files
    elif signal.dtype == 'int32':
        nb_bits = 32 # -> 32-bit wav files
    max_nb_bit = float(2 ** (nb_bits - 1))
    signal = signal / (max_nb_bit + 1.0) # scale the signal to [-1.0,1.0]

    if fs_signal != fs :
        signalr = resample_poly(signal, fs, fs_signal)
    else:
        signalr = signal

    return signalr, fs

def ci_processing(xr,M=16,Ts=0.001,pulsei=2,MCL=[],THL=[],IDR=60):
    # takes a 16 kHz input signal (xr)

    msat = 20*np.log10(2**15 - 1)
    fs=16000; 

    compression = False
    if MCL and THL:
        if len(MCL) == M and len(THL) == M:
            compression = True
        else:
            print('The MCL and THL arrays must have the same size as the value of M - No compression will be applied')

    # read the wavfile and resample the signal to 16kHz
    #xr = wavfile_read(wavfile,fs)

    #create the biphasic pulse
    fc=2/Ts
    Ns=int(Ts*fs)
    ts=np.linspace(0,Ts,num=Ns)
    p1=0.5*square(2*np.pi*fc*ts)
    p2=0.5*square(2*np.pi*fc/2*ts+np.pi)
    pulse=p1+p2
    if pulsei == 1:
        pulse = np.abs(pulse)

    #filterbank
    if M > 1:
        subsignal = np.zeros((M,xr.shape[0]))
        cf = np.logspace(2.5441,3.7404,num=M) # center frequencies - 350 to 5500 Hz
        freql = np.zeros(M)
        freqh = np.zeros(M)
        Nf = 6
        
        # compute the bandpass filter boundaries based on the center frequencies
        for cbi in range(0,M):
            if cbi == 0: # compute boundaries for first filter
                freql[cbi] = cf[0]**2/np.sqrt(cf[0]*cf[1])
                freqh[cbi] = np.sqrt(cf[0]*cf[1])
            elif cbi == M-1: # compute boundaries for last filter
                freql[cbi] = np.sqrt(cf[-2]*cf[-1])
                freqh[cbi] = cf[-1]**2/np.sqrt(cf[-2]*cf[-1])
            else:
                freql[cbi] = np.sqrt(cf[cbi-1]*cf[cbi])
                freqh[cbi] = np.sqrt(cf[cbi]*cf[cbi+1])
            b,a=butter(Nf,[freql[cbi]/(fs/2),freqh[cbi]/(fs/2)],'bandpass') # compute boundaries based on center frequency
            subsignal[cbi,:] = filtfilt(b,a,xr)
            Nft = Nf
            while abs(max(subsignal[cbi,:].min(), subsignal[cbi,:].max(), key=abs)) > abs(max(xr.min(), xr.max(), key=abs)):
                Nf -= 1 # sometimes the filter order needs to be lower to avoid overshoot
                b,a=butter(Nf,[freql[cbi]/(fs/2),freqh[cbi]/(fs/2)],'bandpass') # compute boundaries based on center frequency
                subsignal[cbi,:] = filtfilt(b,a,xr)
            Nf = Nft
    else:
        subsignal=xr
        subsignal.shape = (1,xr.shape[0])

    # compression
    db_ref = 94 #dB reference for converting to dB
    subsignalr = np.zeros((subsignal.shape))
    subsignalp = subsignalr
    Xfilt = np.zeros((M,subsignal.shape[1]))
    Y = Xfilt
    for cbi in range(0,M):
        #half wave rectification - the negative parts of the signal are omitted / the abs is estimated
        #subsignalr[cbi,subsignal[cbi,:]>0]=subsignal[cbi,subsignal[cbi,:]>0]
        subsignalr[cbi,:]=np.abs(subsignal[cbi,:])
        for i in range(0,subsignal.shape[1],Ns):
            Xfilt[cbi,i]=20*np.log10(np.mean(subsignalr[cbi,i:i+Ns]))+db_ref
            if not compression: 
            #no dynamic compression applied
                Y[cbi,i]=Xfilt[cbi,i]
            else: 
            #apply HiRes compression formula
                Y[cbi,i]=(MCL[cbi]-THL[cbi])*(Xfilt[cbi,i]-msat+12+IDR)/IDR + THL[cbi]
            subsignalp[cbi,i:i+Ns] = 10**((Y[cbi,i]-db_ref)/20) *pulse

    #with open('input.txt','wb') as f:
    #    np.savetxt(f, subsignal, fmt='%.3e',delimiter=' ',newline='\n')
    #with open('output.txt','wb') as f:
    #    np.savetxt(f, subsignalp, fmt='%.3e',delimiter=' ',newline='\n')

    return subsignalp

if __name__ == '__main__':
    args = build_argparser().parse_args()

    wavfile = args.wavfile
    M = args.bands
    Ts = args.pulse_duration
    pulsei = args.pulse
    MCL = args.comfortable_levels
    THL = args.threshold
    idr = args.dynamic_range

    fs=16000
    xr = wavfile_read(wavfile,fs)
    subsignalp = ci_processing(xr,M,Ts,pulsei,MCL,THL,IDR)

