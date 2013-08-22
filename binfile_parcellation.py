#### binfile_parcellation.py
# Copyright (C) 2010 R. Cameron Craddock (cameron.craddock@gmail.com)
#
# This script is a part of the pyClusterROI python toolbox for the spatially
# constrained clustering of fMRI data. It performs subject level normalized
# clustering of connectivity matrices. 
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
# python_ncut_lib distributed with this script to be installed in a directory
# that is accessible through PythonPath (the current directory will do for
# python_ncult_lib.py).
import time as time
from numpy import *
from scipy.sparse import csc_matrix
from python_ncut_lib import *

# binfile_parcellate( infile, outfile, K ):
#
# This function performs normalized cut clustering on the connectivity matrix
# specified by infile into sets of K clusters.
#    infile:  .NPY or .bin file containing a representation of the connectivity
#             matrix to be clustered. This file contains a single vector of
#             length 3*N, in which the first N values correspond to the i
#             indices, the second N values correspond to the j indices and the
#             last N values correspond to the weights w_ij of the similarity
#             matrix W. For more information on constructing the input files
#             refer to make_local_connectivity_tcorr.py,
#             make_local_connectivity_scorr.py or
#             make_local_connectivity_ones.py.
#    outfile: a prefix for the output file, this name will be suffixed by
#             _K.npy where K corresponds to the clustering level
#    K:       list of numbers of clusters that will be generated. If this is a
#             single number then only that clustering will be generated. If
#             this is a list of numbers, then the normalized cut algorithm will
#             be run several times, once for each k in the list, and a seperate
#             output file will be generated for each clustering
# 
def binfile_parcellate( infile, outfile, K ):

    # check how long it takes
    start=time.time()
    print 'started at ',start

    # read in the file, I used to use .bin files, but now I use .npy as they
    # contain informaiton about datatype, still support both
    if( infile.endswith(".npy") ):
        print "Reading",infile,"as a npy filetype"
        a=load(infile)
    else:
        print "Reading",infile,"as a binary file of doubles"
        fileobj=open(infile, 'rb')
        a=fromfile(fileobj)
        fileobj.close()

    # calculate the number of non-zero weights in the connectivity matrix
    n=len(a)/3

    # reshape the 1D vector read in from infile in to a 3xN array
    a=reshape(a,(3,n))
    m=max(max(a[0,:]),max(a[1,:]))+1

    # make the sparse matrix, CSC format is supposedly efficient for matrix
    # arithmetic
    W=csc_matrix((a[2,:],(a[0,:],a[1,:])), shape=(m,m))

    print 'finished reading in data and calculating connectivity after ',\
        time.time()-start,'\n'

    # we only have to calculate the eigendecomposition of the LaPlacian once,
    # for the largest number of clusters provided. This provides a significant
    # speedup, without any difference to the results.
    Kmax=max(K)
    eigenval,eigenvec = ncut(W,Kmax)

    print 'finished calculating eigenvectors ',time.time()-start,'\n'

    # calculate each desired clustering result
    for k in K:
        eigk=eigenvec[:,:k]
        eigenvec_discrete = discretisation(eigk)
        print 'finished discretisation ',k,' at ',time.time()-start,'\n'

        # transform the discretised eigenvectors into a single vector
        # where the value corresponds to the cluster # of the corresponding
        # ROI
        group_img=eigenvec_discrete[:,0]

        for i in range(1,k):
            group_img=group_img+(i+1)*eigenvec_discrete[:,i]

        # apply the suffix to the output filename and write out results
        # as a .npy file
        outname=outfile+'_'+str(k)+'.npy'
        save(outname,group_img.todense())

    print 'finished after ',time.time()-start,'\n'
