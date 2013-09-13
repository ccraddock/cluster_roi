#### make_local_connectivity_scorr.py
# Copyright (C) 2010 R. Cameron Craddock (cameron.craddock@gmail.com)
#
# This script is a part of the pyClusterROI python toolbox for the spatially
# constrained clustering of fMRI data. It constructs a spatially constrained
# Coefficient between the whole brain FC maps generaged by the voxel
# timecourses.
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
# NiBabel (http://nipy.sourceforge.net/nibabel) to be installed in a
# directory that is accessible through PythonPath 
import nibabel as nb
from numpy import array
from scipy import *
from scipy.sparse import *
from numpy import corrcoef

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
    if( rank(idx) == 1):
        idx1=idx[0]*prod(sz[1:3])+idx[1]*sz[2]+idx[2]
    else:
        idx1=idx[:,0]*prod(sz[1:3])+idx[:,1]*sz[2]+idx[:,2]
    return idx1

# make_local_connectivity_scorr( infile, maskfile, outfile, thresh )
#
# This script is a part of the ClusterROI python toolbox for the spatially
# constrained clustering of fMRI data. It constructs a spatially constrained
# connectivity matrix from a fMRI dataset. The weights w_ij of the connectivity
# matrix W correspond to the _spatial_correlation_ between the whole brain FC
# maps generated from the time series from voxel i and voxel j. Connectivity is
# only calculated between a voxel and the 27 voxels in its 3D neighborhood
# (face touching and edge touching). The resulting datafiles are suitable as
# inputs to the function binfile_parcellate.
#
#     infile:   name of a 4D NIFTI file containing fMRI data
#     maskfile: name of a 3D NIFTI file containing a mask, which restricts the
#               voxels used in the analysis
#     outfile:  name of the output file, which will be a .NPY file containing
#               a single 3*N vector. The first N values are the i index, the
#               second N values are the j index, and the last N values are the
#               w_ij, connectivity weights between voxel i and voxel j.
#     thresh:   Threshold value, correlation coefficients lower than this value
#               will be removed from the matrix (set to zero).
#
def make_local_connectivity_scorr( infile, maskfile, outfile, thresh ):
    # index
    neighbors=array([[-1,-1,-1],[0,-1,-1],[1,-1,-1],
                     [-1, 0,-1],[0, 0,-1],[1, 0,-1],
                     [-1, 1,-1],[0, 1,-1],[1, 1,-1],
                     [-1,-1, 0],[0,-1, 0],[1,-1, 0],
                     [-1, 0, 0],[0, 0, 0],[1, 0, 0],
                     [-1, 1, 0],[0, 1, 0],[1, 1, 0],
                     [-1,-1, 1],[0,-1, 1],[1,-1, 1],
                     [-1, 0, 1],[0, 0, 1],[1, 0, 1],
                     [-1, 1, 1],[0, 1, 1],[1, 1, 1]])


    # read in the mask
    msk=nb.load(maskfile)
    msz=msk.shape

    # convert the 3D mask array into a 1D vector
    mskdat=reshape(msk.get_data(),prod(msz))

    # determine the 1D coordinates of the non-zero 
    # elements of the mask
    iv=nonzero(mskdat)[0]
    
    # read in the fmri data
    # NOTE the format of x,y,z axes and time dimension after reading
    # nb.load('x.nii.gz').shape -> (x,y,z,t)
    nim=nb.load(infile)
    sz=nim.shape
    print sz, ' dimensions of the 4D fMRI data'
    

    # reshape fmri data to a num_voxels x num_timepoints array
    imdat=reshape(nim.get_data(),(prod(sz[:3]),sz[3]))
    
    # mask the datset to only then in-mask voxels
    imdat=imdat[iv,:]
    imdat_sz = imdat.shape
    
    #zscore fmri time courses, this makes calculation of the
    # correlation coefficient a simple matrix product
    imdat_s=tile(std(imdat,1),(imdat_sz[1],1)).T
    # replace 0 with really large number to avoid div by zero
    imdat_s[imdat_s==0]=1000000
    imdat_m=tile(mean(imdat,1),(imdat_sz[1],1)).T
    imdat=(imdat-imdat_m)/imdat_s
    
    # set values with no variance to zero
    imdat[imdat_s==0]=0
    imdat[isnan(imdat)]=0
    
    # remove voxels with zero variance, do this here
    # so that the mapping will be consistent across
    # subjects
    vndx=nonzero(var(imdat,1)!=0)[0]
    iv=iv[vndx]
    
    m = len(iv)
    print m , ' # of non-zero valued or non-zero variance voxels in the mask'

    # construct a sparse matrix from the mask
    msk=csc_matrix((vndx+1,(iv,zeros(m))),shape=(prod(msz),1))

    sparse_i=[]
    sparse_j=[]
    sparse_w=[[]]

    for i in range(0,m):
        if i % 1000 == 0: print 'voxel #', i
        # convert index into 3D and calculate neighbors
        ndx3d=indx_1dto3d(iv[i],sz[:-1])+neighbors
        # convert resulting 3D indices into 1D
        ndx1d=indx_3dto1d(ndx3d,sz[:-1])
        # convert 1D indices into masked versions
        ondx1d=msk[ndx1d].todense()
        # exclude indices not in the mask
        ndx1d=ndx1d[nonzero(ondx1d)[0]]
        ndx1d=ndx1d.flatten()
        ondx1d=array(ondx1d[nonzero(ondx1d)[0]])
        ondx1d=ondx1d.flatten()-1
        # keep track of the index corresponding to the "seed"
        nndx=nonzero(ndx1d==iv[i])[0]
        # extract the time courses corresponding to the "seed"
        # and 3D neighborhood voxels
        tc=imdat[ondx1d,:]
        # calculate functional connectivity maps for "seed"
        # and 3D neighborhood voxels               
        fc=dot(tc,imdat.T)/(sz[3]-1)
        # calculate the spatial correlation between FC maps
        R=corrcoef(fc)
        if rank(R) == 0:
            R=reshape(R,(1,1))
        # set NaN values to 0
        R[isnan(R)]=0
        # set values below thresh to 0
        R[R<thresh]=0
        # keep track of the indices and the correlation weights
        # to construct sparse connectivity matrix
        sparse_i=append(sparse_i,ondx1d,0)
        sparse_j=append(sparse_j,(ondx1d[nndx])*ones(len(ondx1d)))
        sparse_w=append(sparse_w,R[nndx,:],1)


    # insure that the weight vector is the correct shape
    sparse_w=reshape(sparse_w,prod(shape(sparse_w)))

    # concatenate the i, j, and w_ij vectors
    outlist=sparse_i
    outlist=append(outlist,sparse_j)
    outlist=append(outlist,sparse_w)

    # save the output file to a .NPY file
    save(outfile,outlist)

    print 'finished ',infile,' len ',len(outlist)
