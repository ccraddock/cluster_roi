#### pyClusterROI_test.py
# Copyright (C) 2010 R. Cameron Craddock (cameron.craddock@gmail.com)
#
# This script is a part of the pyClusterROI python toolbox for the spatially
# constrained clustering of fMRI data. It is a demonstration of how to use the
# toolbox and a regression test to make sure that the toolbox code works.
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
# NiBabel (http://niftilib.sourceforge.net/pynifti/) and the pyClusterROI
# toolbox to be installed in a directory that is accessible through PythonPath

# this is how you would add a directory to the search path, this is useful if
# you are running this script from a directory other than the directory where
# the pyClusterROI is installed. Or if for some reason your NumPy, SciPy, or
# NiBabel libraries are in a non-standard location, do this before you import
# the files/libraries that require the change in path
# import sys
# sys.path.append("/home/user/python_toolboxes")

# import the different functions we will use from pyClusterROI

# only need one of these, based on which connectivity metric you prefer
from make_local_connectivity_ones import *
from make_local_connectivity_scorr import *
from make_local_connectivity_tcorr import *

# do not need this if you are peforming group mean clustering
from binfile_parcellation import *

# import the functions for group clustering, only need one of these
from group_binfile_parcellation import *
from group_mean_binfile_parcellation import *

# import if you want to write the results out to nifti, only need
# one of these, probably just want the one that does renumbering,
# why do i include the other one? no idea.
from make_image_from_bin import *
from make_image_from_bin_renum import *

from time import time

T0 = time()

# the name of the maskfile that we will be using
maskname="gm_maskfile.nii.gz"

# make a list of all of the input fMRI files that we will be using
infiles = [  'subject1.nii.gz', 'subject2.nii.gz', 'subject3.nii.gz' ]

##### Step 1. Individual Conenctivity Matrices 
# first we need to make the individual connectivity matrices, I will
# do this for all three different kinds (tcorr, scorr, ones) but you
# will only need to do it for one

# the easiest is random clustering which doesn't require any functional
# data, just the mask
print 'ones connectivity'
make_local_connectivity_ones( maskname, 'rm_ones_connectivity.npy')


# construct the connectivity matrices using tcorr and a r>0.5 threshold

for idx, in_file in enumerate(infiles):

    # construct an output filename for this file
    outname='rm_tcorr_conn_'+str(idx)+'.npy'

    print 'tcorr connectivity',in_file
    # call the funtion to make connectivity
    make_local_connectivity_tcorr( in_file, maskname, outname, 0.5 )


# construct the connectivity matrices using scorr and a r>0.5 threshold
# This can take a _really_ long time
for idx, in_file in enumerate(infiles):

    # construct an output filename for this file
    outname='rm_scorr_conn_'+str(idx)+'.npy'

    print 'scorr connectivity',in_file
    # call the funtion to make connectivity
    make_local_connectivity_scorr( in_file, maskname, outname, 0.5 )

##### Step 2. Individual level clustering
# next we will do the individual level clustering, this is not performed for 
# group-mean clustering, remember that for these functions the output name
# is a prefix that will have K and .npy added to it by the functions. We
# will perform this for clustering between 100, 150 and 200 clusters
NUM_CLUSTERS = [100,150,200]

# For random custering, this is all we need to do, there is no need for group
# level clustering, remember that the output filename is a prefix, and 
binfile_parcellate('rm_ones_connectivity.npy','rm_ones_cluster',NUM_CLUSTERS)

# for tcorr
for idx, in_file in enumerate(infiles):

    # construct filenames
    infile='rm_tcorr_conn_'+str(idx)+'.npy'
    outfile='rm_tcorr_indiv_cluster_'+str(idx)

    print 'tcorr parcellate',in_file
    binfile_parcellate(infile, outfile, NUM_CLUSTERS)

# for scorr
for idx, in_file in enumerate(infiles):

    # construct filenames
    infile='rm_scorr_conn_'+str(idx)+'.npy'
    outfile='rm_scorr_indiv_cluster_'+str(idx)

    print 'scorr parcellate',in_file
    binfile_parcellate(infile, outfile, NUM_CLUSTERS)

##### Step 3. Group level clustering
# perform the group level clustering for clustering results containing 100, 150,
# and 200 clusters. as previously mentioned, this does _not_ have to be done for
# random clustering

# for both group-mean and 2-level clustering we need to know the number of
# nonzero voxels in in the mask 
mask_voxels=(nb.load(maskname).get_data().flatten()>0).sum()

# group_mean clustering is pretty simple, input the connectivity files and run.
# we can perform multiple clusterings with this function, so once again the
# output filename is a prefix
tcorr_conn_files=['rm_tcorr_conn_0.npy','rm_tcorr_conn_1.npy',\
    'rm_tcorr_conn_2.npy']
print 'group-mean parcellate tcorr'
group_mean_binfile_parcellate( tcorr_conn_files,\
    'rm_group_mean_tcorr_cluster', NUM_CLUSTERS,mask_voxels);

# now group mean cluster scorr files
scorr_conn_files=['rm_scorr_conn_0.npy','rm_scorr_conn_1.npy',\
    'rm_scorr_conn_2.npy']
print 'group-mean parcellate scorr'
group_mean_binfile_parcellate( scorr_conn_files,\
    'rm_group_mean_scorr_cluster', NUM_CLUSTERS, mask_voxels);

# the 2-level clustering has to be performed once for each desired clustering
# level, and requires individual level clusterings as inputs
for k in NUM_CLUSTERS:
    ind_clust_files=[]
    for i in range(0,len(infiles)):
        ind_clust_files.append('rm_tcorr_indiv_cluster_'+str(i)+\
            '_'+str(k)+'.npy')

    print '2-level parcellate tcorr',k
    group_binfile_parcellate(ind_clust_files,\
        'rm_group_tcorr_cluster_'+str(k)+'.npy',k,mask_voxels)

# now for scorr 
for k in NUM_CLUSTERS:
    ind_clust_files=[]
    for i in range(0,len(infiles)):
        ind_clust_files.append('rm_scorr_indiv_cluster_'+str(i)+\
            '_'+str(k)+'.npy')

    print '2-level parcellate scorr',k
    group_binfile_parcellate(ind_clust_files,\
        'rm_group_scorr_cluster_'+str(k)+'.npy',k,mask_voxels)

##### Step 4. Convert the binary output .npy files to nifti
# this can be done with or without renumbering the clusters to make sure they
# are contiguous. remember, we might end up with fewer clusters than we ask for,
# and this could result in gaps in the cluster numbering. Choose which you like,
# i use them intermittently below as a regression test

# write out for the random clustering

for k in NUM_CLUSTERS:
    binfile='rm_ones_cluster_'+str(k)+'.npy'
    imgfile='rm_ones_cluster_'+str(k)+'.nii.gz'
    make_image_from_bin(imgfile,binfile,maskname);

# write out for group mean clustering
for k in NUM_CLUSTERS:
    binfile='rm_group_mean_tcorr_cluster_'+str(k)+'.npy'
    imgfile='rm_group_mean_tcorr_cluster_'+str(k)+'.nii.gz'
    make_image_from_bin_renum(imgfile,binfile,maskname)

for k in NUM_CLUSTERS:
    binfile='rm_group_mean_scorr_cluster_'+str(k)+'.npy'
    imgfile='rm_group_mean_scorr_cluster_'+str(k)+'.nii.gz'
    make_image_from_bin_renum(imgfile,binfile,maskname)

# write out for group 2-level clustering
for k in NUM_CLUSTERS:
    binfile='rm_group_tcorr_cluster_'+str(k)+'.npy'
    imgfile='rm_group_tcorr_cluster_'+str(k)+'.nii.gz'
    make_image_from_bin_renum(imgfile,binfile,maskname)

for k in NUM_CLUSTERS:
    binfile='rm_group_scorr_cluster_'+str(k)+'.npy'
    imgfile='rm_group_scorr_cluster_'+str(k)+'.nii.gz'
    make_image_from_bin_renum(imgfile,binfile,maskname)

T1 = time()

print '******************************'
print 'time is ', T1-T0
##### FIN

