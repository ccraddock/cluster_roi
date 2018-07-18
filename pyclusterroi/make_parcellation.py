# make_parcellation.py
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
import nibabel as nb
import numpy as np
import time as time
import python_ncut_lib as nc
from .make_connectivity import *
import scipy.sparse as sp
import datetime


def discretize_and_write(eigenvectors, nim, mask_array, out_file):

    # prepare the output nifti file

    if len(eigenvectors.shape) == 2:
        num_clusters = eigenvectors.shape[1]
    else:
        raise ValueError("Expected eigenvectors to be a 2D matrix, received {0}-D".format(len(eigenvectors.shape)))

    discrete_eigenvectors = nc.discretisation(eigenvectors)

    # transform the discrete eigenvectors into a single vector
    # where the value corresponds to the cluster # of the corresponding
    # ROI
    cluster_image = discrete_eigenvectors[:, 0].todense()

    for k in range(1, num_clusters):
        cluster_image += (k + 1) * discrete_eigenvectors[:, k].todense()

    # remove any gaps in numbering
    cluster_image_nogap = cluster_image.copy()
    i = 1
    for j in np.unique(cluster_image.tolist()):
        cluster_image_nogap[cluster_image == j] = i
        i += 1

    # setup array for outputting images
    image_data = np.zeros((np.prod(nim.shape)), dtype='uint16')

    # store the output
    image_data.flat[mask_array > 0] = np.short(cluster_image_nogap.flatten())

    thdr = nim.header
    thdr['scl_slope'] = 1
    nim_aff = nim.affine

    # reshape the image
    nim_out = nb.Nifti1Image(image_data.reshape(nim.shape), nim_aff, thdr)

    nim_out.set_data_dtype('uint16')
    nim_out.to_filename(out_file)

    return out_file


def make_parcellation(num_clusters, method, image_files, mask_file, out_file, thresh=0.5, num_threads=1):

    """
    Make connectivity for an individual fMRI scan contained in a nifti image file, or the mean connectivity matrices
    calculated across several fMRI scans. This can be used for a variety of scenarios including both individual and
    group level clustering.

    (1) To obtain individual level clustering from a single or multiple fMRI datasets (e.g. averaging multiple datasets
      from the same individual), provide a list of one or more 4D fMRI datasets and choose tcorr or scorr for the method

    (2) To obtain group level clustering by averaging connectivity graphs across groups and then parcellating the
      result, provide a list of one or more 4D fMRI datasets and choose tcorr or scorr for the method

    (3) To perform a 2-level clustering, where each individual dataset is parcellated and the results are fed into a
      group level parcellation, first perform (1) for each of the input datasets and then call this method again
      with a list of the 1st level parcellation results to be included, and the 'cluster' method

    (4) To perform a clustering that only uses information about geometry, omit the imaging files and set
      method = "cluster". The parcellation will use the geometry information in mask_file.

    If more than one image_file are provided, the result will be the average of the resulting adjacency matrices. In
    this case the outputs will always be weighted even if the connectivity method (described below) results in binary
    matrices.

    After calculating connectivity, performs ncut parcellation procedure on an adjacency matrix, this procedure
    involves:
    1. normalizing the adjacency matrix
    2. calculating the top eigenvectors of the normalized adjacency matrix, to save resources we
       only calculate the eigenvectors for the largest number of clusters that are requested, solutions for lower
       number of clusters will use the same (albeit fewer) eigenvectors

    and for each cluster solution requested:
    3. use a discretization procedure to rotate eigenvectors into cluster assignments vectors
    4. write out the results

    These last two can be done in parallel for differing cluster solutions.


    :param num_clusters: number of desired clusters to be returned, if a list, a file will be written for a clustering
        result calculated for each num_cluster (k)
    :param method: string indicating how connectivity should be calculated
              should be one of 'cluster', 'tcorr', 'scorr', or 'ones', which are defined as:

              cluster: construct a binary adjacency matrix from a parcellation or atlas where two voxels are considered
                 connected if they are both in the same cluster, and not connected otherwise
              ones: construct a binary adjacency matrix from a mask file where two voxels are considered connected if
                 they are touching neighbors (face and edge)
              tcorr: construct a binary adjacency matrix from 4D fMRI data, where two voxels are considered connected
                 if they are touching neighbors (face and edge) and the pearson's correlation between their time series
                 are greater than thresh
              scorr: construct a binary adjacency matrix from 4D fMRI data, where two voxels are considered connected
                 if they are touching neighbors (face and edge) and the pearson's correlation between their whole brain
                 functional connectivity maps (calculated from the temporal correlation between a voxel's time series
                 and the time series of every voxel in the brain) are greater than thresh

    :param image_files: path to a single input nifti file or a list of paths to several input files, used to calculate
            connectivity, if  a list, connectivity matrices are averaged across files, if method == 'ones' this can
            safely be set to None
    :param mask_file: mask for confining the neuroimaging data to an area of interest, most
              commonly the grey matter
    :param out_file: filename for the output file, occurrences of the '{k}' substring will be replaced by the number
            of clusters in the parcellation, otherwise this will be appended to the out_file string before the
            .nii.gz
    :param thresh: thresh hold to apply to connectivity measures
    :param num_threads: The last stage of clustering involves mapping individual voxels to their cluster assignments.
           This is done once for each of the clustering solutions requested. Although this isn't the most consumptive
           part of the process, it is the easiest to parallelize when the user requests multiple cluster solutions.
           This parameter determines the max number of parallel processes may be used for this process. The actual
           level of parallelization that is used will be determined by the smaller of this number and the number of
           clustering solutions requested.

    :return: 0 on success, raises exceptions on errors
    """

    if isinstance(image_files, str):
        image_files = [image_files]

    start_time = time.time()

    print('Started reading files and calculating connectivity at {0}'.format(
        datetime.datetime.fromtimestamp(start_time).strftime('%H:%M:%S')))

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

    if method in ['tcorr', 'scorr', 'cluster']:

        if image_files:

            adjacency_matrix = None
            for image_file in image_files:

                # load and conform the fMRI data
                image_array = nb.load(image_file).get_data()

                num_volumes = 0

                if 'cluster' in method:

                    num_volumes = 1

                    if len(image_array.shape) > 3 and image_array.shape[3] > 1:
                        print("Warning: Cluster connectivity requires only one volume, received {0}, "
                              "just using the first".format(image_array.shape[1]))

                        image_array = image_array[:, 0]

                    elif len(image_array.shape) < 3 or len(image_array.shape) > 4:
                        raise ValueError(
                            "Cluster connectivity expects a 3 or 4 -D file, received {0}, which is {1}-D".format(
                                image_file, len(image_array.shape)))

                elif method in ['tcorr', 'scorr']:

                    if len(image_array.shape) == 4 and image_array.shape[3] > 3:
                        num_volumes = image_array.shape[3]

                    elif len(image_array.shape) != 4:
                        raise ValueError("Expecting 4-D file {0} is {1}-D.".format(
                            image_file, len(image_array.shape)))

                    elif image_array.shape[3] < 4:
                        raise ValueError("{0} connectivity requires > 3 volumes, {1} has {2}".format(
                            method, image_file, image_array.shape[3]))

                # reshape fmri data to a num_voxels x num_volumes array
                image_array = np.reshape(image_array, (np.prod(image_array.shape[:3]), num_volumes))

                # reduce fmri data to just in-brain voxels
                image_array = image_array[mask_array != 0, :]

                if 'cluster' in method:
                    (w, i, j) = make_local_connectivity_clusters(image_array)
                else:
                    (w, i, j) = make_local_connectivity(method, image_array, np.reshape(mask_array, mask_shape), thresh,
                                                        num_threads)

                # make the sparse matrix, CSC format is supposedly efficient for matrix
                # arithmetic
                if isinstance(adjacency_matrix, type(None)):
                    adjacency_matrix = 1.0 / len(image_files) * \
                               sp.csc_matrix((w, (i, j)), shape=(num_mask_voxels, num_mask_voxels))
                else:
                    adjacency_matrix += 1.0 / len(image_files) * \
                                sp.csc_matrix((w, (i, j)), shape=(num_mask_voxels, num_mask_voxels))

        else:
            raise ValueError(
                "Calculating local connectivity using {0} requires at least one input dataset, got NONE".format(method))

    elif 'ones' in method:
        (w, i, j) = make_local_connectivity(method, None, np.reshape(mask_array, mask_shape), thresh, num_threads)

        # make the sparse matrix, CSC format is supposedly efficient for matrix
        # arithmetic
        adjacency_matrix = sp.csc_matrix((w, (i, j)), shape=(num_mask_voxels, num_mask_voxels))

    else:
        raise ValueError("unknown connectivity method {0}".format(method))

    print('Finished reading in data and calculating connectivity after {0} seconds'.format(
        time.time() - start_time))

    if isinstance(num_clusters, int):
        num_clusters = [num_clusters]

    start_time = time.time()

    # load mask
    nim = nb.load(mask_file)

    mask_array = ((nim.get_data()) > 0).astype('uint8')
    mask_array = mask_array.flatten()

    # we only have to calculate the eigen decomposition of the LaPlacian once,
    # for the largest number of clusters provided. This provides a significant
    # speedup, without any difference to the results.
    k_max = max(num_clusters)
    eigenvalues, eigenvectors = nc.ncut(adjacency_matrix, k_max)

    print('Finished calculating eigenvectors after {0} seconds'.format(time.time() - start_time))

    # if outfile includes a file extension, remove it
    out_file = out_file.replace(".nii", "")
    out_file = out_file.replace(".gz", "")

    if '{k}' not in out_file:
        out_file = out_file + '_variant-{method}{k:04d}_roi'

    out_file += ".nii.gz"

    out_files = []

    if num_threads == 1:

        for k in num_clusters:
            out_files.append(discretize_and_write(eigenvectors[:, :k], nim, mask_array,
                                                  out_file.format(**{'method':method, 'k': k})))

    elif num_threads > 1:
        import multiprocessing as mp

        # calculate each desired clustering result
        with mp.Pool(processes=num_threads) as pool:
            apply_results = [
                pool.apply_async(discretize_and_write,
                                 (eigenvectors[:, :k], nim, mask_array, out_file.format(**{'k': k})))
                for k in num_clusters]

            out_files = [result.get() for result in apply_results]

    else:
        raise ValueError("Expected num_threads to be an integer greater than zero, received {0}.".format(num_threads))

    return out_files
