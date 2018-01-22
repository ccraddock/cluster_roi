#### make_parcellation.py
# Copyright (C) 2010 R. Cameron Craddock (cameron.craddock@gmail.com)
#
# This script is a part of the pyClusterROI python toolbox for the spatially
# constrained clustering of fMRI data. It computes ncut clustering on
# connectivity matrices calculated using the functions in make_connectivity.py
#
# For more information refer to:
#
# Craddock, R. C.; James, G. A.; Holtzheimer, P. E.; Hu, X. P. & Mayberg, H. S.
# A whole brain fMRI atlas generated via spatially constrained spectral
# clustering Human Brain Mapping, 2012, 33, 1914-1928 doi: 10.1002/hbm.21333.
#
# ARTICLE{Craddock2012,
#   author = {Craddock, R C and James, G A and Holtzheimer, P E and Hu, X P and
#   Mayberg, H S},
#   title = {{A whole brain fMRI atlas generated via spatially constrained
#   spectral clustering}},
#   journal = {Human Brain Mapping},
#   year = {2012},
#   volume = {33},
#   pages = {1914--1928},
#   number = {8},
#   address = {Department of Neuroscience, Baylor College of Medicine, Houston,
#       TX, United States},
#   pmid = {21769991},
# }
#
# Documentation, updated source code and other information can be found at the
# NITRC web page: http://www.nitrc.org/projects/cluster_roi/ and on github at
# https://github.com/ccraddock/cluster_roi
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
####
from .make_connectivity import *


def make_individual_parcellation(num_clusters, method, image_file, mask_file, out_file, thresh=0.5):
    """
    Make connectivity for an individual fMRI scan contained in a nifti image file"

    :param num_clusters: number of desired clusters to be returned, can be a list
    :param method: string indicating how connectivity should be calculated
              should be one of tcorr, scorr, or ones
    :param image_file: nifti file containing neuroimaging data to calculate connectivity
              from
    :param mask_file: mask for confining the neuroimaging data to an area of interest, most
              commonly the grey matter
    :param out_file: output nifti file containing the parcellated data
    :param thresh: thresh hold to apply to connectivity measures

    :return: no return, raises exceptions on errors
    """

    import nibabel as nb
    import scipy.sparse as sp
    import numpy as np
    import time as time
    import python_ncut_lib as nc

    if isinstance(num_clusters, int):
        num_clusters = [num_clusters]

    start_time = time.time()

    print 'started at ', start_time

    # load mask
    nim = nb.load(mask_file)
    mask_array = ((nim.get_data()) > 0).astype('uint32')

    mask_shape = mask_array.shape
    mask_array = mask_array.flatten()

    # make the mask do joint duty as a lookup table
    num_mask_voxels = 0
    for mask_index, mask_value in enumerate(mask_array):
        if mask_value == 1:
            num_mask_voxels += 1
            mask_array[mask_index] = num_mask_voxels

    # load and conform the fMRI data
    image_array = nb.load(image_file).get_data()

    image_shape = image_array.shape
    num_time_courses = image_shape[-1]

    # reshape fmri data to a num_voxels x num_time_points array
    image_array = np.reshape(image_array, (np.prod(image_shape[:-1]), num_time_courses))

    # reduce fmri data to just in-brain voxels
    image_array = image_array[mask_array != 0, :]

    print("parcellating im_array: {0} {1}".format(image_array.shape[0], image_array.shape[1]))

    if 'tcorr' in method:
        (w, i, j) = make_local_connectivity_tcorr(image_array, np.reshape(mask_array, mask_shape), thresh)

    elif 'scorr' in method:
        (w, i, j) = make_local_connectivity_scorr(image_array, np.reshape(mask_array, mask_shape), thresh)

    elif 'ones' in method:
        (w, i, j) = make_local_connectivity_ones(np.reshape(mask_array, mask_shape))

    else:
        raise ValueError("unknown connectivity method {0}".format(method))

    # make the sparse matrix, CSC format is supposedly efficient for matrix
    # arithmetic
    w_matrix = sp.csc_matrix((w, (i, j)), shape=(num_mask_voxels, num_mask_voxels))

    print('finished reading in data and calculating connectivity after {0} seconds'.format(
        time.time() - start_time))

    # we only have to calculate the eigendecomposition of the LaPlacian once,
    # for the largest number of clusters provided. This provides a significant
    # speedup, without any difference to the results.
    k_max = max(num_clusters)
    eigenvalues, eigenvectors = nc.ncut(w_matrix, k_max)

    print('finished calculating eigenvectors after ({0}) seconds'.format(time.time() - start_time))

    # setup array for outputting images
    image_data = np.zeros((np.prod(mask_shape)), dtype='uint16')

    # if outfile includes a file extension, remove it
    out_file = out_file.replace(".nii", "")
    out_file = out_file.replace(".gz", "")

    # prepare the output nifti file
    thdr = nim.get_header()
    thdr['scl_slope'] = 1

    nim_aff = nim.get_affine()

    # calculate each desired clustering result
    for k in num_clusters:
        discrete_eigenvectors = nc.discretisation(eigenvectors[:, :k])
        print('Finished discretisation {0} after {1} seconds.'.format(k, time.time()-start_time))

        # transform the discrete eigenvectors into a single vector
        # where the value corresponds to the cluster # of the corresponding
        # ROI
        cluster_image = discrete_eigenvectors[:, 0].todense()

        for i in range(1, k):
            cluster_image += (i+1) * discrete_eigenvectors[:, i].todense()

        # remove any gaps in numbering
        cluster_image_nogap = cluster_image.copy()
        i = 1
        for j in np.unique(cluster_image.tolist()):
            cluster_image_nogap[cluster_image == j] = i
            i += 1

        # store the output
        image_data[mask_array > 0] = np.short(cluster_image_nogap.flatten())

        # reshape the image
        nim_out = nb.Nifti1Image(image_data.reshape(mask_shape), nim_aff, thdr)

        nim_out.set_data_dtype('uint16')
        nim_out.to_filename("{0}_nclust-{1:04d}.nii.gz".format(out_file, k))

    return out_file