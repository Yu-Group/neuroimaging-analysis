import numpy as np

from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import ball, dilation
from skimage.measure import label, regionprops
from skimage.io import imread

from sklearn.neighbors import KDTree

def apply_threshold(im, thresh=None, min_thresh=0, gauss_kernel_sz=1, max_quant=.999, verbose=False):
    '''Filter out noise below some intensity threshold.

    Params
    ------
    im: ndarray, a 3D volume
    thresh: if None (default), uses Otsu's method to calculate intensity threshold
    min_thresh: lower threshold of voxel intensity
    gauss_kernel_size: 3D Guassian kernel size to smooth image
    max_quant: if thresh is None, use as max quantile of voxel intensities
    verbose: if True, print the calculated threshold from Otsu's method
    '''
    im_max = np.quantile(im, max_quant)
    if im_max <= 0:
        raise ValueError(f'The {max_quant} quantile of input image intensity is <= 0!')
    im_gauss = im.copy()
    if thresh is None:
        # calculate a threshold using Otsu's method
        # I believe 'mirror' is similar to 'symmetric' as used in ABC's code
        # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter.html
        im_gauss = gaussian(im_gauss, gauss_kernel_sz, mode='mirror', preserve_range=True)
        thresh = threshold_otsu(im[np.logical_and(im > 0, im < im_max)])
        if verbose:
            print(f'Otsu\'s threshold: {thresh:.3f}')
        if thresh < min_thresh:
            thresh = min_thresh
    im_gauss = im_gauss - thresh
    im_gau
    ss[im_gauss < 0] = 0
    return im_gauss

def filter_small_artifacts(im, min_vox_vol=516):
    '''Label the connected components of the image, and filter out cc's with fewer
    pixels than min_vox_vol.

    NOTE: this function has side effects -- it modifies the input image!

    Params
    ------
    im: ndarray, a 3D volume
    min_vox_vol: the minimum number of active voxels an artifact can have

    '''
    labels, num_ccs = label(im > 0, return_num=True)
    cc_sizes = np.array([np.sum(labels == i) for i in range(1, num_ccs+1)])
    invalid_labels = np.arange(1, num_ccs+1)[cc_sizes <= min_vox_vol]
    im[np.isin(labels, invalid_labels)] = 0
    return im

def local_max_3d(im, wdims, clear_border=False):
    '''Find local maxima in 3D.
    Adapted from XR_locmax3d from https://github.com/abcucberkeley/LLSM3DTools.

    Params
    ------
    im: ndarray, a 3D volume
    wdims: list [x, y, z], the dimensions of the local neighborhood to search
    clear_border: if True, set borders to zero
    '''
    if len(wdims) == 1:
        wx = wy = wz = wdims[0]
    elif len(wdims) == 3:
        wx, wy, wz = wdims[0], wdims[1], wdims[2]
    else:
        raise ValueError('wdims should be a list with either 1 or 3 numbers')

    im_cp = im.copy()
    im_cp[np.isnan(im)] = 0
    dilated_im = dilation(im_cp, np.ones((wx, wy, wz)))

    # if input is equal to max filter response, then the point is a local max
    im_cp[dilated_im != im_cp] = 0

    # set borders to zero
    if clear_border:
        hwx = wx // 2
        hwy = wy // 2
        hwz = wz // 2
        im_cp[:(hwx - 1), :, :] = 0
        im_cp[-hwx:, :, :] = 0
        im_cp[:, :(hwy - 1), :] = 0
        im_cp[:, -hwy:, :] = 0
        im_cp[:, :, :(hwz - 1)] = 0
        im_cp[:, :, -hwz:] = 0

    return im_cp

def remove_duplicates(pts, min_dist):
    '''Remove points that are within min_dist of one another, replacing them
    with their mean.

    NOTE: this function has side effects -- it modifies the input array!

    Params
    ------
    pts: 2d ndarray, points to filter
    min_dist: the minimum distance to consider points the same
    '''
    tree = KDTree(pts)
    dists, inds = tree.query(pts, k=2)
    dup_rows = dists[:, 1] <= min_dist

    if np.sum(dup_rows) == 0:
        return pts

    pts_cp = pts.copy()

    # I don't think the following line (reproduced from MATLAB) is necessary
    # dup_inds = tree.query(pts_data[dup_rows, :], k=2, return_distance=False)

    dup_pts = pts_cp[dup_rows, :]
    dup_inds = inds[dup_rows, :]
    dup_dists = dists[dup_rows, 1]

    delete_idx = []

    while np.sum(dup_rows) > 0:
        idx = dup_inds[np.argmin(dup_dists), 0] # the index of point closest to another
        nn_idx = inds[idx, 1] # the nearest neighbors index
        pts_cp[idx] = np.around((pts_cp[idx, :] + pts_cp[nn_idx, :]) / 2)

        # indices of dup_pts whose nearest neighbors are one of these two pts
        dup_idx = np.logical_and(~np.isin(dup_inds[:, 0], (idx, nn_idx)),
                                 np.isin(dup_inds[:, 1], (idx, nn_idx)))

        if not np.all(~dup_idx):
            # reset the distances given the new point
            dists[dup_inds[dup_idx, 0]] = [np.linalg.norm(pts_cp[idx] - pt) for pt in dup_pts[dup_idx]]
            inds[dup_inds[dup_idx, 0], 1] = idx

        # remove the point at nn_idx
        delete_idx.append(nn_idx)
        dists[idx] = np.Inf
        dists[nn_idx] = np.Inf

        dup_rows = dists[:, 1] <= min_dist
        dup_pts = pts_cp[dup_rows, :]
        dup_inds = inds[dup_rows, :]
        dup_dists = dists[dup_rows, 1]

    # recurse in case some new points fall within the min_dist of each other
    return remove_duplicates(np.delete(pts_cp, delete_idx, axis=0), min_dist)

def calc_synapse_density(im,
                         expansion_factor = 4.09,
                         px_sz = 0.097,
                         px_sz_z = 0.18,
                         sigmas = [2.5, 2.5],
                         density_vol = 1,
                         find_local_max = True,
                         min_az = 0.1,
                         thresh = None,
                         min_thresh = 0,
                         gauss_kernel_sz = 1,
                         max_quant = .999,
                         min_vox_vol = 516,
                         ball_radius = 2,
                         min_nonzero_frac = 0.1,
                         verbose = False):
    '''Calculate the synapse density in a given volume. Adapted from Gokul
    Upadhyayula's MATLAB code for the same task.

    Params
    ------
    im: ndarray, a 3D volume
    expansion_factor: expansion factor by which the vol expanded from original size
    px_sz: pixel size in x,y plane, in micrometers
    px_sz_z: pixel size in z direction, in micrometers
    sigmas: for 3D gaussian filter
    density_vol: pre-expanded surveyed volume, in micrometers cubed
    find_local_max: whether to use min_az to remove duplicates (TODO: of what?)
    min_az: threshold in nanometers (pre-expanded) for active zone duplicates,
            only used when find_local_max is True
    thresh: if None (default), uses Otsu's method to calculate intensity threshold
    min_thresh: lower threshold of voxel intensity
    gauss_kernel_sz: 3D Guassian kernel size to smooth image
    max_quant: max quantile of data for Otsu's method calculation
    min_vox_vol: set to zero if no "pre-cleaning" is necessary to remove
                 high-freq noise
    ball_radius: to make a ball structure element from the local max
    min_nonzero_frac: the loaded volume must have this fraction of voxels with values > 0

    Returns
    -------
    A dictionary with the following keys:
      synapse_density_mask, clean_vol, synapse_indices, synapse_densities
    '''
    bw = min_az / 2 / px_sz * expansion_factor # search radius
    rad = np.cbrt(density_vol * 3 / 4 / np.pi)
    rad_ex = rad * expansion_factor
    bw_ex = rad_ex / px_sz # search radius after expansion
    struc_el = ball(ball_radius)
    px_z_ratio = px_sz_z / px_sz

    if np.sum(im) > 0 and np.count_nonzero(im) >= (im.size * min_nonzero_frac):
        im_gauss = apply_threshold(im, thresh, min_thresh, gauss_kernel_sz, max_quant)
    else:
        ValueError(f'There weren\'t enough non-zero voxels in the input volume!')

    if min_vox_vol > 0:
        # filter out small artifacts from the image
        im_gauss = filter_small_artifacts(im_gauss)

    if np.sum(im_gauss > 0) == 0:
        ValueError(f'There weren\'t any non-zero voxels after thresholding!')

    # calculate local maxima or centroids of synapses
    if find_local_max:
        all_max = local_max_3d(
            im_gauss, (2*np.ceil([sigmas[0], sigmas[0], sigmas[1]])+1).astype(int),
            clear_border=False
        )
        max_ind = np.nonzero(all_max)
        puncta_ind = np.array([max_ind[0], max_ind[1], max_ind[2]]).T
    else:
        centroids = regionprops_table(labels, properties=['centroid'])
        puncta_ind = np.array([centroids['centroid-0'],
                               centroids['centroid-1'],
                               centroids['centroid-2']]).T

    if len(puncta_ind[0]) == 1:
        ValueError(f'There was only one synapse centroid/maximum!')

    if find_local_max:
        # filter detections closer than 2*bw, replacing with the mean of close points
        puncta_ind = remove_duplicates(puncta_ind, 2*bw)

    # stretch out the z axis
    puncta_ind_z = puncta_ind.copy()
    puncta_ind_z[:, 2] = puncta_ind_z[:, 2] * px_z_ratio

    # get the number of puncta_ind within a radius of bw_ex of each point and
    # divide by density_vol to calculate the synapse density
    tree = KDTree(puncta_ind_z)
    neighbor_count = tree.query_radius(puncta_ind_z, r=bw_ex, count_only=True)
    synapse_density = neighbor_count / density_vol

    mask = np.zeros(im.shape)
    mask[puncta_ind[:, 0], puncta_ind[:, 1], puncta_ind[:, 2]] = synapse_density
    mask = dilation(mask, struc_el)

    nonzero_idx = im_gauss != 0
    clean_vol = np.zeros(im.shape)
    clean_vol[nonzero_idx] = im[nonzero_idx]

    out_dict = dict()
    out_dict['synapse_density_mask'] = mask
    out_dict['clean_vol'] = clean_vol
    out_dict['synapse_ind'] = synapse_ind
    out_dict['synapse_density'] = synapse_density

    return out_dict
