#### make_local_connectivity_ones.py
# Copyright (C) 2010 R. Cameron Craddock (cameron.craddock@gmail.com)
#
# This script is a part of the pyClusterROI python toolbox for the spatially
# constrained clustering of fMRI data. It constructs a spatially constrained
# connectivity matrix from a fMRI dataset, where then connectivity weight
# betwen two neighboring voxels is 1. This implements "random" clustering
# described in the paper below.
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

# this scripts requires NumPy (numpy.scipy.org), SciPy (www.scipy.org), and
# NiBabel (http://nipy.sourceforge.net/nibabel) to be installed in a directory
# that is accessible through PythonPath 
import nibabel as nb
import numpy as np
#from numpy import array
from scipy import *
from scipy.sparse import *

# simple function to translate 1D vector coordinates to 3D matrix coordinates,
# for a 3D matrix of size sz
def indx_1dto3d(idx,sz):
    x=divide(idx,prod(sz[1:3]))
    y=divide(idx-x*prod(sz[1:3]),sz[2])
    z=idx-x*prod(sz[1:3])-y*sz[2]
    return (x,y,z)

# simple function to translate 3D matrix coordinates to 1D vector coordinates,
# for a 3D matrix of size sz
def indx_3dto1d(idx,sz):
    if( ndim(idx) == 1):
        idx1=idx[0]*prod(sz[1:3])+idx[1]*sz[2]+idx[2]
    else:
        idx1=idx[:,0]*prod(sz[1:3])+idx[:,1]*sz[2]+idx[:,2]
    return idx1

# array to add to 3d coordinates of seed voxel to get the face, edge, and corner
# touching neighbors
neighbors = np.array([[-1, -1, -1], [0, -1, -1], [1, -1, -1],
                      [-1, 0, -1], [0, 0, -1], [1, 0, -1],
                      [-1, 1, -1], [0, 1, -1], [1, 1, -1],
                      [-1, -1, 0], [0, -1, 0], [1, -1, 0],
                      [-1, 0, 0], [0, 0, 0], [1, 0, 0],
                      [-1, 1, 0], [0, 1, 0], [1, 1, 0],
                      [-1, -1, 1], [0, -1, 1], [1, -1, 1],
                      [-1, 0, 1], [0, 0, 1], [1, 0, 1],
                      [-1, 1, 1], [0, 1, 1], [1, 1, 1]])

# make_local_connectivity_ones( maskfile, outfile )
#
# This script is a part of the ClusterROI python toolbox for the spatially
# constrained clustering of fMRI data. It constructs a spatially constrained
# connectivity matrix for a fMRI dataset. The weights w_ij of the connectivity
# matrix W are set to 1 if a voxel is withen the 3D neighborhood (face touching
# and edge touching) of the center voxel. The resulting datafiles are suitable
# as inputs to the function binfile_parcellate.
#
#     maskfile: name of a 3D NIFTI file containing a mask, which restricts the
#               voxels used in the analysis
#     outfile:  name of the output file, which will be a .NPY file containing
#               a single 3*N vector. The first N values are the i index, the
#               second N values are the j index, and the last N values are the
#               w_ij, connectivity weights between voxel i and voxel j.
#
def make_local_connectivity_ones( maskfile ):

    # read in the mask
    msk=nb.load(maskfile)
    msz=msk.shape
    mskdat=msk.get_data().flatten()
    ix=np.nonzero(mskdat)[0]
    print '%d non-zero voxels in the mask'%(len(ix))

    # get the range of the ix's
    ix_min = np.min(ix)
    ix_max = np.max(ix)

    mvals=[]

    # loop over all of the voxels in the mask 	
    for i in ix:
        print i
    	if i % 1000 == 0: print 'voxel #', i

        seed_neigh = [x for x in indx_3dto1d(indx_1dto3d(i, msz) + neighbors, msz) if
                      ix_min <= x <= ix_max and mskdat[x] != 0]
        mvals+=zip(int(i)*np.ones(len(seed_neigh),dtype=np.int),seed_neigh,np.ones(len(seed_neigh)))


    return(mvals)
