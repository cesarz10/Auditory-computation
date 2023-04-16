"""
Signal Preprocessing functions.
Sarah Verhulst, Deepak Baby, UGent, March 2019.
"""

import numpy as np
import math

def apply_preemphasis (x, coef):
    '''
    Applies preemphasis filtering on signal x 
    coef is the preemphasis filter coefficient
    '''
    return np.append(x[0], x[1:] - coef * x[:-1])

def frame_timedomain (x, framelength, frameshift):
    '''
    Split a time-domain signal x into frames of length framelength
    shifed by frameshift
    Both framelength and frameshift must be specified as sample lengths
    '''
    sig_length = len(x)
    num_frames = ((sig_length - framelength) // frameshift) + 1
    indices = np.tile(np.arange(framelength),(num_frames, 1)) + np.tile(np.arange(0, num_frames * frameshift, frameshift), (framelength, 1)).T 
    frames = x[indices.astype(np.int32, copy=False)]
    return frames

def apply_hamming(x):
    '''
    Apply hamming window on the framed signal x
    '''
    winlength = x.shape[1]
    return x * np.hamming(winlength)

def pow_stft(x, framelength, frameshift, nfft):
    '''
    Computes power short-time Fourier spectrogram of a signal x
    With framelength and frameshift with nfft frequency bins
    '''
    frames = frame_timedomain(x, framelength, frameshift) # split signal into frames
    frames = apply_hamming(frames) # apply hamming window
    mag_stft = np.absolute(np.fft.rfft(frames, nfft))
    return (mag_stft ** 2) / nfft

def adjust_spl(x, level):
    '''
    Adjust the soud pressure level of a signal x to the desired level
    '''
    p0 = 2e-5 # reference SPL
    target_energy = p0 * (10 ** (level/20)) # scaling required for level SPL
    x_rms = np.sqrt(np.mean(np.square(x))) # rms of the signal
    return x * target_energy / x_rms


def next_power_of_two (n):
    '''
    Find the smallest power of 2 greater than n
    '''
    return 1 if n == 0 else 2**math.ceil(math.log(n) / math.log(2))

def cochleagram(x, winlength, winshift):
    """
    Converts a time-domain filtered signal into energies
    expects a 2D input x
    """
    nchan, x_len = x.shape
    n_frames  = ((x_len - winlength) // winshift) + 1
    out = np.zeros((nchan, n_frames))
    ham = np.hamming(winlength)
    for m in range(n_frames):
        startpoint = m * winshift
        endpoint = startpoint + winlength
        sig_segment = np.multiply(x[:, startpoint : endpoint], ham)
        out[:,m] = np.sum(np.square(sig_segment), axis=1)
    return out

def lifter_mfcc(mfcc):
    cep_lifter = 23 # popularly used value
    (nframes, ncoeff) = mfcc.shape
    n = np.arange(ncoeff)
    lift = 1 + (cep_lifter / 2) * np.sin(np.pi * n / cep_lifter)
    return mfcc * lift