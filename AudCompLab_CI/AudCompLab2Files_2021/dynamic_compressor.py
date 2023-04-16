#!/usr/bin/env python

"""
Dynamic Range Compressor implemented from 
"Implementation and evaluation of an experimental hearing aid dynamic
range compressor" by Giso Grimm et al.
 
Fotis Drakopoulos, UGent, March 2020

signal, fsy: input signal to be processed and its sampling frequency
(default = 16 kHz)- the signal will be downsampled to 16 kHz for the
processing.
 
Gtable(k,L) is the gain table fitted for an individual hearing loss 
profile. Each element of the table corresponds to an insertion gain 
applied for each respective frequency and for each respective input 
level. The values can be either positive (amplification) or negative 
(attenuation). The 1st dimension of the table corresponds to the 
frequency bins of the STFT of the input (given in the freq(k) array) and 
the 2nd dimension corresponds to the input levels in dB-SPL (given in the
gt(L) array). If no argument is given, zero gain is applied.
The Gtable can also be a column vector. In this case, the gain applied is
only frequency dependent (the same gain is applied for all input levels).
The Gtable always needs to be two-dimensional, therefore the second 
dimension (level) will be 1 in this case. In the same way, the Gtable can 
be a row vector for a frequency-independent gain prescription.

freq: The frequencies of the STFT bins corresponding to the 1st dimension
of the Gtable. If not defined, the default frequencies for each bin are 
computed by f=(1:(W))*(16e3/(2*W)), where W = 64 by default.
gt: The input intensity levels in dB-SPL corresponding to the 2nd
dimension of the Gtable. The values need to be positive.
 
W: frame size (default = 64) - twice the size is used for the FFT (after 
zero-padding) and 50% overlap is used between the frames
ta: attack time constant (default = 20ms)
tr: release time constant (default = 100ms)

spl_norm variable is used to define the scaling of the magnitudes to
dB-SPL (default is ~94)
"""

from __future__ import division, print_function
import sys
from argparse import ArgumentParser
import numpy as np
from scipy.signal import square, butter, filtfilt, stft, istft, resample_poly
from scipy.interpolate import interp1d
import scipy.io.wavfile

def build_argparser():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_signal", help="Input signal to be processed", required=True, type=float)
    parser.add_argument("-fs", "--fs_signal", help="Sampling frequency of the input signal", required=True, type=int)
    parser.add_argument("-G", "--Gain_table", help="Gain prescription table to be applied", default=[], type=float)
    parser.add_argument("-gt", "--level_axis", help="Level axis of the gain table in dB-SPL", default=[], type=float)
    parser.add_argument("-f", "--freq_axis", help="Frequency axis of the gain table in Hz", default=[], type=float)
    parser.add_argument("-W", "--frame_size", help="Frame size of the analysis", default=64, type=int)
    parser.add_argument("-ta", "--attack", help="Attack time constant in secs", default=0.02, type=float)
    parser.add_argument("-tr", "--release", help="Release time constant in secs", default=0.1, type=float)

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

def dynamic_compression(signal,fs_signal=16000,Gtable=[],freq=[],gt=[],W=64,ta=0.02,tr=0.1):

    # check the inputs
    if isinstance(Gtable,list):
        Gtable=np.array(Gtable)
    if Gtable.size != 0 and Gtable.ndim != 2:
        print("The gain table must be always two-dimensional (frequency x level)")
        sys.exit()
    if Gtable.ndim == 2 and Gtable.size > 0:
        if Gtable.shape[0] != freq.shape[0] or Gtable.shape[1] != gt.shape[0]:
            print("The dimensions of the gt and freq don't match the dimensions of the Gtable provided")
            sys.exit()

    if gt.size < 2 and Gtable.shape[1]==1:
        gtstep=10. # gain table granulation in dB
        gtmin=10. # minimum is 10 dB-SPL
        gtmax=100. # maximum is 100 dB-SPL
        gt=np.arange(gtmin,gtmax,gtstep)
        if Gtable.size==0:
            Gtable = np.zeros((W,gt.size));
            print('No gain table provided - No compression will be applied')
        else:
            Gtable = np.tile(Gtable,(1,gt.size))

    spl_norm=20*np.log10(2e-5) # dB reference of the spectrogram to convert the scale to dB-SPL
    db_offset=200 # just an offset to perform logarithmic interpolation over negative values

    h=0.5 # overlap percentage 
    npad=W # zero-padding size 
    fs=16e3 # processing sampling frequency
    if fs_signal != fs :
        xr = resample_poly(signal, fs, fs_signal) # resample the signal
    else:
        xr = signal

    pads=int((1-h)*W)
    sp=int(W*h)
    '''
    wnd=hamming(W)
    x=np.pad(xr,(pads,pads))

    # segment signal for FFT (STFT)
    L=x.shape[0]
    N=int((L-W)/sp+1) # number of segments
    Index=(np.tile(np.arange(1,W),(N,1))+np.tile(np.arange(0,(N-1))*sp,(1,W)))
    hw=np.tile(wnd,(1,N))
    xs=x[Index]*hw

    # pad with zeros
    x_pad=np.pad(xs,((npad/2,npad/2),(0,0)))

    # STFT
    X=np.fft.fft(xs) # whole STFT

    XPhase=np.angle(X[1:np.ceil(X.shape[0]/2),:]) # Phase
    X2=np.abs(X[1:np.ceil(X.shape[0]/2),:])**2; # Spectrogram
    '''
    f, _, X=stft(xr,fs=fs,window='hann',nperseg=W,noverlap=sp,nfft=W+2*pads,return_onesided=True)
    XPhase=np.angle(X) # Phase
    X2=(np.abs(X))**2 # Spectrogram
    frames=X.shape[1]

    if freq.size < 2 and Gtable.shape[0]==1:
        freq=f
        Gtable = np.tile(Gtable,(freq.size,1))
    else:
        # if the gain table provided is not described for each frequency bin,
        # compute the rest of them via logarithmic interpolation
        Gtablet=Gtable
        Gtable=np.zeros((f.shape[0],Gtablet.shape[1]))
        for n in range (0,gt.shape[0]):
            interp_fun = interp1d(np.log10(freq),Gtablet[:,n],kind='linear',fill_value='extrapolate')
            Gtable[:,n] = interp_fun(np.log10(f))
        freq=f
        del Gtablet, interp_fun
    del f
    ind=np.isnan(Gtable)
    Gtable[ind]=0
    ind=np.isinf(Gtable)
    Gtable[ind]=0

    aa=np.exp(-1./(ta*fs)) # attack coefficient
    ar=np.exp(-1./(tr*fs)) # release coefficient

    G=np.zeros(X2.shape)
    Lin=np.zeros(X2.shape)
    Xlev=np.zeros(X2.shape)
    Lin[:,0]=10*np.log10(X2[:,0])-spl_norm # input intensity - first bin
    for n in range(1,frames):
        # input intensity for the gain table - positive and negative bins are
        # taken into account
        Xlev[:,n]=10*np.log10(2*X2[:,n])-spl_norm
        for k in range(0,Xlev.shape[0]):
            # apply the attack and release filters to the input intensity
            if Xlev[k,n] <= Lin[k,n-1]:
                Lin[k,n] = aa*Lin[k,n-1]+(1-aa)*Xlev[k,n]
            else:
                Lin[k,n] = ar*Lin[k,n-1]+(1-ar)*Xlev[k,n]

            if Lin[k,n]<0: # negative values are not needed
                Lin[k,n]=0

            interp_fun = interp1d(np.log10(gt+db_offset),Gtable[k,:],kind='linear',fill_value='extrapolate')
            G[k,n] = interp_fun(np.log10(Lin[k,n]+db_offset))
            # gain to be applied computed with logarithmic interpolation
            # extrapolation used for values outside the given SPL range.

    Y2=10*np.log10(X2)+G # processed STFT

    #y=OverlapAdd2(10.^(Y2/20),XPhase,(W+npad),h*W);
    #y=y(pads+npad/2+1:end-pads-npad/2); # processing signal in time domain

    Y=(10**(Y2/20))*np.exp(1j*XPhase)
    _, y=istft(Y,fs=fs,window='hann',nperseg=W,noverlap=sp,nfft=W+2*pads,input_onesided=True)

    if fs_signal != fs :
        signalp = resample_poly(y, fs_signal, fs) # resample the signal to the original fs
    else:
        signalp = y

    return signalp

if __name__ == '__main__':
    args = build_argparser().parse_args()

    signal=args.input_signal
    fs_signal=args.fs_signal
    Gtable=args.Gain_table
    gt=args.level_axis
    freq=args.freq_axis
    W=args.frame_size
    ta=args.attack
    tr=args.release

    signalp = dynamic_compression(signal,fs_signal,Gtable,gt,freq,W,ta,tr)

