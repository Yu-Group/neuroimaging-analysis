import re
import time
import pickle
import numpy as np

from tqdm import tqdm
from os import listdir
from os.path import join as oj
from cv2 import imread, IMREAD_GRAYSCALE

from matplotlib import gridspec
import matplotlib.pyplot as plt

def call_and_timing(fn, *args, _fn_name: str=None, _verbose: bool=True, **kwargs):
    start = time.time()
    result = fn(*args, **kwargs)
    duration = time.time() - start
    if _verbose:
        if _fn_name is None:
            print(f'Duration of call: {duration:.3}s')
        else:
            print(f'Duration of "{_fn_name}" call: {duration:.3}s')
    return result, duration

def read_vol(path_to_slices, n_slices=None, first_slice=0, verbose=False):
    if n_slices is None:
        regex = re.compile(r'\d+') 
        n_slices = np.amax([int(regex.findall(f)[0]) for f in listdir(path_to_slices) if '.tif' in f]) + 1
    
    if verbose:
        print(f'Reading {n_slices} slices from {path_to_slices}, starting at {first_slice}.tif')
        seq = tqdm(range(1, n_slices))
    else:
        seq = range(1, n_slices)
    
    # read first slice to get shape
    slice0 = imread(oj(path_to_slices, f'{first_slice}.tif'), IMREAD_GRAYSCALE)
    
    # initialize data structure for all slices
    vol = np.empty((n_slices,) + slice0.shape)
    
    # add first slice
    vol[0, :, :] = slice0
    
    # add remaining slices        
    for i in seq:
        vol[i, :, :] = imread(oj(path_to_slices, f'{first_slice+i}.tif'), IMREAD_GRAYSCALE)
        
    if verbose:
        print(f'Volume shape: {vol.shape}')

    return vol
    
def save_as_pickle(obj, filepath):
    with open(filepath, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def plot_channels(chs: list, figsize=None, max_quant=.999):
    R = 1
    C = len(chs)
    if C > 3:
        R = C // 3
        C = 3
    figsize = (C+1, R+1) if figsize is None else figsize
    plt.subplots(figsize=figsize, sharey=True)
    gs = gridspec.GridSpec(R, C,
             wspace=0.0, hspace=0.0, 
             top=1.-0.5/(R+1), bottom=0.5/(R+1), 
             left=0.5/(C+1), right=1-0.5/(C+1))
    
    for r in range(R):
        for c in range(C):
            ch_cp = chs[r*3 + c].copy()
            ch_max = np.quantile(ch_cp, 0.999)
            ch_cp[ch_cp > ch_max] = ch_max
            ax = plt.subplot(gs[r,c])
            ax.imshow(ch_cp, vmin=0, vmax=ch_max)

    plt.show()
