import re
import time
import pickle
import numpy as np

from tqdm import tqdm
from os import listdir
from os.path import join as oj
from cv2 import imread, IMREAD_GRAYSCALE

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
