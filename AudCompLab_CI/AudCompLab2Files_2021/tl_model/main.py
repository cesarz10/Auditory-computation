"""
Main function to load the data and do the various analyses.
Sarah Verhulst, Deepak Baby, UGent, 2019
"""

import numpy as np
import scipy.io as sio
from tl_model import tl_cochlea

# read the wavfile
fs, x = sio.wavfile.read('example.wav')