### Import packages
import time as time
import numpy as np
from numpy import *
from scipy import *
from scipy.sparse import *
import nibabel as nb
import os
import sys
from collections import defaultdict
from scipy.linalg import norm, svd, LinAlgError

### nifti_gen_fix
"""
 Temporary band-aid script to generate NIfTI images from NumPY parcellation data.
 Written by Dan Lurie (danjlurie@gmail.com).
 Usage:
 python nifti_gen_fix.py /path/to/parcel_data.npy /path/to/mask.nii.gz /path/to/new_atlas_image.nii.gz
 """
def nifti_gen_fix(parcel_path, mask_path, out_path):
    # Read in command line arguments.
    parcel_path, mask_path, out_path = sys.argv[1:4]
    # Load parcellation data.
    parcel_data = np.load(parcel_path)
    # Load the mask file.
    mask_file = nb.load(mask_path)
    # Get the mask data.
    mask_data = mask_file.get_data()
    # Recast the mask data as float 64 to match parcellation data.
    atlas_data = mask_data.astype('float64')
    # Copy the parcellation data to the within-mask voxels.
    atlas_data[atlas_data == 1] = parcel_data.ravel()
    # Get the affine and header from the mask.
    mask_affine = mask_file.get_affine()
    mask_header = mask_file.get_header()
    # Set the data type for the header to float64.
    mask_header.set_data_dtype('float64')
    # Create a NIfTI image from the atlas data.
    atlas_image = nib.Nifti1Image(atlas_data, mask_affine, mask_header)
    # Write the NIfTI file
    atlas_image.to_filename(out_path)
    print("Fixed")

### Simple functions
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

# exception hander for singular value decomposition
class SVDError(Exception):
    def __init__(self,value):
        self.value=value
    def __str__(self):
        return repr(self.value)

### binfile_parcellation
"""
binfile_parcellate( infile, outfile, K ):

This function performs normalized cut clustering on the connectivity matrix
specified by infile into sets of K clusters.
   infile:  .NPY or .bin file containing a representation of the connectivity
            matrix to be clustered. This file contains a single vector of
            length 3*N, in which the first N values correspond to the i
            indices, the second N values correspond to the j indices and the
            last N values correspond to the weights w_ij of the similarity
            matrix W. For more information on constructing the input files
            refer to make_local_connectivity_tcorr.py,
            make_local_connectivity_scorr.py or
            make_local_connectivity_ones.py.
   outfile: a prefix for the output file, this name will be suffixed by
            _K.npy where K corresponds to the clustering level
   K:       list of numbers of clusters that will be generated. If this is a
            single number then only that clustering will be generated. If
            this is a list of numbers, then the normalized cut algorithm will
            be run several times, once for each k in the list, and a seperate
            output file will be generated for each clustering
"""
def binfile_parcellate( infile, outfile, K ):
    # check how long it takes
    start=time.time()
    print('started at ',start)
    # read in the file, I used to use .bin files, but now I use .npy as they
    # contain informaiton about datatype, still support both
    if( infile.endswith(".npy") ):
        print("Reading",infile,"as a npy filetype")
        a=load(infile)
    else:
        print("Reading",infile,"as a binary file of doubles")
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
    print('finished reading in data and calculating connectivity after ',\
        time.time()-start,'\n')
    # we only have to calculate the eigendecomposition of the LaPlacian once,
    # for the largest number of clusters provided. This provides a significant
    # speedup, without any difference to the results.
    Kmax=max(K)
    eigenval,eigenvec = ncut(W,Kmax)
    print('finished calculating eigenvectors ',time.time()-start,'\n')
    # calculate each desired clustering result
    for k in K:
        eigk=eigenvec[:,:k]
        eigenvec_discrete = discretisation(eigk)
        print('finished discretisation ',k,' at ',time.time()-start,'\n')
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
    print('finished after ',time.time()-start,'\n')

### group_binfile_parcellation
"""
group_img=group_binfile_parcellate( infiles, outfile, K,n_voxels ):

This function performs group level clustering of individual level clustering
results. Each single subject clustering is converted into a coincidence matrix
W, where w_ij = 1 if voxels i and j are in the same cluster, and zero
otherwise. Coincidence matrices are averaged across subjects and then
submitted to normalized cut clustering to obtain K clusters.
   infiles:   list of .NPY or .bin file containing single subject clustering
              results. Each of these files is a 1D vector where the ith value
              corresponds to the cluster assignment of voxel i. This assumes
              that the vectors are in the same order. These files are
              generated by binfile_parcellation.py
   outfile:   the name of the output file, a n_voxels x 1 vector is written to
              this file, where the ith value corresponds to the cluster to
              which voxel i is assigned
   K:         The number of clusters that will be generated. This is a single
              number. This assumes that each of the input files were clustered
              to K
   n_voxels:  Number of voxels in the _mask_ used to generate the subject
              specific connectivity matrices
   group_img: (output) n_voxels x 1 vector indicating cluster assignments
"""
def group_binfile_parcellate( infiles, outfile, K, n_voxels ):
    # read in all of the files, construct coincidence matrices, and average
    # them
    for i in range(0,len(infiles)):
        # read in the files
        if infiles[i].endswith(".npy"):
            print("Reading",infiles[i],"as a npy filetype")
            imdat=load(infiles[i])
        else:
            print("Reading %s as a binary file of doubles"%(\
                infiles[i]))
            fid=open(infiles[i], 'rb')
            imdat=fromfile(fid)
            fid.close()
        # construct the coincidences between the voxels
        sparse_i=[]
        sparse_j=[]
        for j in range(1,K+1):
            grp_ndx=nonzero(imdat==j)[0]
            ff=tile(grp_ndx,(len(grp_ndx),1))
            sparse_i=append(sparse_i,reshape(ff,prod(shape(ff))))
            sparse_j=append(sparse_j,reshape(ff.transpose(),prod(shape(ff))))
        # sum the coincidence matrices across input files
        if i==0:
            W=csc_matrix((ones(len(sparse_i)),(sparse_i,sparse_j)),\
               (n_voxels,n_voxels),dtype=double)
        else:
            print('adding ',i)
            W=W+csc_matrix((ones(len(sparse_i)),(sparse_i,sparse_j)),\
                (n_voxels,n_voxels),dtype=double)
    # set the diagonal to zeros as is customary with affinity matrices
    W=W-spdiags(W.diagonal(),[0],n_voxels,n_voxels,"csc")
    print("diag is ",sum(W.diagonal()),"\n")
    # divide by the number of input files, to calculate the average
    W=W/len(infiles)
    print("finished reading in data and calculating connectivity\n")
    # perform clustering
    eigenval,eigenvec = ncut(W,K)
    eigenvec_discrete = discretisation(eigenvec)
    ## write out the results
    group_img=eigenvec_discrete[:,0]
    for i in range(1,K):
    	if not i%10: print(i)
        group_img=group_img+(i+1)*eigenvec_discrete[:,i]
    save(outfile,group_img.todense())
    print("finished group parcellation\n")
    return array(group_img.todense())

### group_mean_binfile_parcellation
"""
group_mean_binfile_parcellate( infiles, outfile, K ):

This function performs group level clustering of individual level clustering
results. Each single subject clustering is converted into a coincidence
matrix W, where w_ij = 1 if voxels i and j are in the same cluster, and zero
otherwise. Coincidence matrices are averaged across subjects and then
submitted to normalized cut clustering to obtain K clusters.
   infiles:   list of .NPY or .bin file containing single subject clustering
              results. Each of these contains a 3*N x 1 vector where the
              first N values coorespond to the i indices, the next N
              correspond to the j indices and the last N values correspond to
              the weights w_ij. For more information on constructing these
              files refer to make_local_connectivity_tcorr.py,
              make_local_connectivity_scorr.py, or
              make_local_connectivity_ones.py.
   outfile:   a prefix for the output file, this name will be suffixed by
              _K.npy where K corresponds to the clustering level
   K:         list of numbers of clusters that will be generated. If this is
              a single number then only that clustering will be generated. If
              this is a list of numbers, then the normalized cut algorithm
              will be run several times, once for each k in the list, and a
              seperate output file will be generated for each clustering
   n_voxels:  Number of voxels in the _mask_ used to generate the subject
              specific connectivity matrices
"""
def group_mean_binfile_parcellate( infiles, outfile, K, n_voxels ):
    if not infiles or not outfile or not K or not n_voxels or K == 0 or n_voxels == 0:
        print("Invalid arguments")
        raise ValueError
    # index
    start=time.time()
    print('started at ',start)
    # read in the files, convert them to similarity matrices,
    # and then average them
    for i in range(0,len(infiles)):
        # read in the file
        if infiles[i].endswith(".npy"):
            print("Reading",infiles[i],"as a npy filetype")
            a=load(infiles[i])
        else:
            print("Reading",infiles[i],"as a binary file of doubles")
            fileobj=open(infiles[i], 'rb')
            a=fromfile(fileobj)
            fileobj.close()
        n=len(a)/3
        a=reshape(a,(3,n))
        # determine all of the voxel indices represented
        vx_ndx=unique(a[-2,:])
        # make the sparse matrix, CSC format is supposedly efficient for matrix
        # arithmetic
        if i==0:
            W=csc_matrix((a[2,:],(a[0,:],a[1,:])), shape=(n_voxels,n_voxels))
        else:
            print('adding ',i)
            W=W+csc_matrix((a[2,:],(a[0,:],a[1,:])), shape=(n_voxels,n_voxels))
    # complete the average
    W=W/len(infiles)
    vx_ndx=unique(W.nonzero()[0])
    print('finished reading in data and calculating connectivity after ',\
        time.time()-start,'\n')
    # we only have to calculate the eigendecomposition of the LaPlacian once,
    # for the largest number of clusters provided. This provides a significant
    # speedup, without any difference to the results.
    Kmax=max(K)
    eigenval,eigenvec = ncut(W,Kmax)
    print('finished calculating eigenvectors ',time.time()-start,'\n')
    # calculate each desired clustering result
    for k in K:
        eigk=eigenvec[:,:k]
        eigenvec_discrete = discretisation(eigk)
        print('finished discretisation ',k,' at ',time.time()-start,'\n')
        # transform the discretised eigenvectors into a single vector
        # where the value corresponds to the cluster # of the corresponding
        # ROI
        group_img=eigenvec_discrete[:,0]
        for i in range(1,k):
            group_img=group_img+(i+1)*eigenvec_discrete[:,i]
        # write out the results
        outname=outfile+'_'+str(k)+'.npy'
        save(outname,group_img.todense())
        print('finished ',k,' after ',time.time()-start,'\n')
    print('finished after ',time.time()-start,'\n')

### make_image_from_bin_renum
"""
make_image_from_bin_renum( image, binfile, mask ):

Converts a NPY file generated by binfile_parcellation.py,
group_binfile_parcellation.py, or group_mean_binfile_parcellation.npy into a
nifti file where each voxels intensity corresponds to the number of the
cluster to which it belongs. Clusters are renumberd to be contiguous.
    image:   The name of the nifti file to be written.
    binfile: The binfile to be converted. The file contains a n_voxel x 1
             vector that is converted to a nifti file.
    mask:    Mask describing the space of the nifti file. This should
             correspond to the mask originally used to create the
             connectivity matrices used for parcellation.
"""
def make_image_from_bin_renum( image, binfile, mask ):
    # read in the mask
    nim=nb.load(mask)
    # read in the binary data
    if( binfile.endswith(".npy") ):
        print("Reading",binfile,"as a npy filetype")
        a=load(binfile)
    else:
        print("Reading",binfile,"as a binary file of doubles")
        a=fromfile(binfile)
    unique_a=list(set(a.flatten()))
    unique_a.sort()
    # renumber clusters to make the contiguous
    b=zeros((len(a),1))
    for i in range(0,len(unique_a)):
        b[a==unique_a[i]]=i+1
    imdat=nim.get_data()
    imdat=imdat.astype('int16')
    # map the binary data to mask
    imdat[imdat>0]=1
    imdat[imdat>0]=short(b[0:sum(imdat)].flatten())
    # Get the mask header and change the dtype
    nim_head = nim.get_header()
    nim_head.set_data_dtype('int16')
    # write out the image as nifti
    nim_out = nb.Nifti1Image(imdat, nim.get_affine(), nim_head)
    #nim_out.set_data_dtype('int16')
    nim_out.to_filename(image)

### make_image_from_bin
"""
make_image_from_bin( image, binfile, mask )

Converts a NPY file generated by binfile_parcellation.py,
group_binfile_parcellation.py, or group_mean_binfile_parcellation.npy into a
nifti file where each voxels intensity corresponds to the number of the
cluster to which it belongs.
    image:   The name of the nifti file to be written.
    binfile: The binfile to be converted. The file contains a n_voxel x 1
             vector that is converted to a nifti file.
    mask:    Mask describing the space of the nifti file. This should
             correspond to the mask originally used to create the
             connectivity matrices used for parcellation.
"""
def make_image_from_bin( image, binfile, mask ):
    # read in the mask
    nim=nb.load(mask)
    # read in the binary data
    if( binfile.endswith(".npy") ):
        print("Reading",binfile,"as a npy filetype")
        a=load(binfile)
    else:
        print("Reading",binfile,"as a binary file of doubles")
        a=fromfile(binfile)
    imdat=nim.get_data()
    print("shape",shape(a))
    print("sum",sum(imdat))
    imdat=imdat.astype('int16')
    # map the binary data to mask
    mask_voxels=(imdat.flatten()>0).sum()
    print("shape2",shape(a[0:mask_voxels]))
    imdat[imdat>0]=short(a[0:mask_voxels].flatten())
    # write out the image as nifti
    thdr=nim.get_header()
    thdr['scl_slope']=1
    thdr.set_data_dtype('int16')
    nim_aff = nim.get_affine()
    nim_out = nb.Nifti1Image(imdat, nim_aff, thdr)
    #nim_out.set_data_dtype('int16')
    nim_out.to_filename(image)

### make_local_connectivity_ones
"""
make_local_connectivity_ones( maskfile, outfile )

This script is a part of the ClusterROI python toolbox for the spatially
constrained clustering of fMRI data. It constructs a spatially constrained
connectivity matrix for a fMRI dataset. The weights w_ij of the connectivity
matrix W are set to 1 if a voxel is withen the 3D neighborhood (face touching
and edge touching) of the center voxel. The resulting datafiles are suitable
as inputs to the function binfile_parcellate.

    maskfile: name of a 3D NIFTI file containing a mask, which restricts the
              voxels used in the analysis
    outfile:  name of the output file, which will be a .NPY file containing
              a single 3*N vector. The first N values are the i index, the
              second N values are the j index, and the last N values are the
              w_ij, connectivity weights between voxel i and voxel j.
"""
def make_local_connectivity_ones( maskfile, outfile ):
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
    m=len(iv)
    print(m, ' # of non-zero voxels in the mask')
    # construct a sparse matrix from the mask
    msk=csc_matrix((range(1,m+1),(iv,zeros(m))),shape=(prod(msz),1))
    sparse_i=[]
    sparse_j=[]
    sparse_w=[]
    # loop over all of the voxels in the mask
    for i in range(0,m):
        if i % 1000 == 0: print('voxel #', i)
        # calculate the voxels that are in the 3D neighborhood
        # of the center voxel
        ndx3d=indx_1dto3d(iv[i],msz)+neighbors
        ndx1d=indx_3dto1d(ndx3d,msz)
        # restrict the neigborhood using the mask
        ondx1d=msk[ndx1d].todense()
        ndx1d=ndx1d[nonzero(ondx1d)[0]]
        ndx1d=ndx1d.flatten()
        ondx1d=array(ondx1d[nonzero(ondx1d)[0]])
        ondx1d=ondx1d.flatten()
        # determine the index of the seed voxel in the neighborhood
        nndx=nonzero(ndx1d==iv[i])[0]
        # the connections between neighbors are all = 1
        R=ones((len(ndx1d),len(ndx1d)))
        if rank(R) == 0:
            R=reshape(R,(1,1))
        # extract just the weights connected to the seed
        R=R[nndx,:].flatten()
        # determine the non-zero correlations (matrix weights)
        # and add their indices and values to the list
        nzndx=nonzero(R)[0]
        if(len(nzndx)>0):
            sparse_i=append(sparse_i,ondx1d[nzndx]-1,0)
            sparse_j=append(sparse_j,(ondx1d[nndx]-1)*ones(len(nzndx)))
            sparse_w=append(sparse_w,R[nzndx],1)
    # concatenate the i, j and w_ij into a single vector
    outlist=sparse_i
    outlist=append(outlist, sparse_j)
    outlist=append(outlist, sparse_w)
    # save the output file to a .NPY file
    save(outfile,outlist)

### make_local_connectivity_scorr
"""
make_local_connectivity_scorr( infile, maskfile, outfile, thresh )

This script is a part of the ClusterROI python toolbox for the spatially
constrained clustering of fMRI data. It constructs a spatially constrained
connectivity matrix from a fMRI dataset. The weights w_ij of the connectivity
matrix W correspond to the _spatial_correlation_ between the whole brain FC
maps generated from the time series from voxel i and voxel j. Connectivity is
only calculated between a voxel and the 27 voxels in its 3D neighborhood
(face touching and edge touching). The resulting datafiles are suitable as
inputs to the function binfile_parcellate.

    infile:   name of a 4D NIFTI file containing fMRI data
    maskfile: name of a 3D NIFTI file containing a mask, which restricts the
              voxels used in the analysis
    outfile:  name of the output file, which will be a .NPY file containing
              a single 3*N vector. The first N values are the i index, the
              second N values are the j index, and the last N values are the
              w_ij, connectivity weights between voxel i and voxel j.
    thresh:   Threshold value, correlation coefficients lower than this value
              will be removed from the matrix (set to zero).
"""
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
    print(sz, ' dimensions of the 4D fMRI data')
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
    print(m , ' # of non-zero valued or non-zero variance voxels in the mask')
    # construct a sparse matrix from the mask
    msk=csc_matrix((vndx+1,(iv,zeros(m))),shape=(prod(msz),1))
    sparse_i=[]
    sparse_j=[]
    sparse_w=[[]]
    for i in range(0,m):
        if i % 1000 == 0: print('voxel #', i)
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
    print('finished ',infile,' len ',len(outlist))

### make_local_connectivity_tcorr
"""
make_local_connectivity_tcorr( infile, maskfile, outfile, thresh )

This script is a part of the ClusterROI python toolbox for the spatially
constrained clustering of fMRI data. It constructs a spatially constrained
connectivity matrix from a fMRI dataset. The weights w_ij of the connectivity
matrix W correspond to the _temporal_correlation_ between the time series
from voxel i and voxel j. Connectivity is only calculated between a voxel and
the 27 voxels in its 3D neighborhood (face touching and edge touching). The
resulting datafiles are suitable as inputs to the function
binfile_parcellate.

    infile:   name of a 4D NIFTI file containing fMRI data
    maskfile: name of a 3D NIFTI file containing a mask, which restricts the
              voxels used in the analysis
    outfile:  name of the output file, which will be a .NPY file containing
              a single 3*N vector. The first N values are the i index, the
              second N values are the j index, and the last N values are the
              w_ij, connectivity weights between voxel i and voxel j.
    thresh:   Threshold value, correlation coefficients lower than this value
              will be removed from the matrix (set to zero).
"""
def make_local_connectivity_tcorr( infile, maskfile, outfile, thresh ):
    # index array used to calculate 3D neigbors
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
    msz=shape(msk.get_data())
    # convert the 3D mask array into a 1D vector
    mskdat=reshape(msk.get_data(),prod(msz))
    # determine the 1D coordinates of the non-zero
    # elements of the mask
    iv=nonzero(mskdat)[0]
    m=len(iv)
    print(m, '# of non-zero voxels in the mask')
    # read in the fmri data
    # NOTE the format of x,y,z axes and time dimension after reading
    # nb.load('x.nii.gz').shape -> (x,y,z,t)
    nim=nb.load(infile)
    sz=nim.shape
    # reshape fmri data to a num_voxels x num_timepoints array
    imdat=reshape(nim.get_data(),(prod(sz[:3]),sz[3]))
    # construct a sparse matrix from the mask
    msk=csc_matrix((range(1,m+1),(iv,zeros(m))),shape=(prod(sz[:-1]),1))
    sparse_i=[]
    sparse_j=[]
    sparse_w=[]
    negcount=0
    # loop over all of the voxels in the mask
    for i in range(0,m):
        if i % 1000 == 0: print('voxel #', i)
        # calculate the voxels that are in the 3D neighborhood
        # of the center voxel
        ndx3d=indx_1dto3d(iv[i],sz[:-1])+neighbors
        ndx1d=indx_3dto1d(ndx3d,sz[:-1])
        #acd
        # restrict the neigborhood using the mask
        ondx1d=msk[ndx1d].todense()
        ndx1d=ndx1d[nonzero(ondx1d)[0]]
        ndx1d=ndx1d.flatten()
        ondx1d=array(ondx1d[nonzero(ondx1d)[0]])
        ondx1d=ondx1d.flatten()
        # determine the index of the seed voxel in the neighborhood
        nndx=nonzero(ndx1d==iv[i])[0]
	#print(nndx,)
	#print(nndx.shape)
	#print(ndx1d.shape)
	#print(ndx1d)
        # exctract the timecourses for all of the voxels in the
        # neighborhood
        tc=matrix(imdat[ndx1d,:])
        # make sure that the "seed" has variance, if not just
        # skip it
        if var(tc[nndx,:]) == 0:
            continue
        # calculate the correlation between all of the voxel TCs
        R=corrcoef(tc)
        if rank(R) == 0:
            R=reshape(R,(1,1))
        # extract just the correlations with the seed TC
        R=R[nndx,:].flatten()
        # set NaN values to 0
        R[isnan(R)]=0
        negcount=negcount+sum(R<0)
        # set values below thresh to 0
        R[R<thresh]=0
        # determine the non-zero correlations (matrix weights)
        # and add their indices and values to the list
        nzndx=nonzero(R)[0]
        if(len(nzndx)>0):
            sparse_i=append(sparse_i,ondx1d[nzndx]-1,0)
            sparse_j=append(sparse_j,(ondx1d[nndx]-1)*ones(len(nzndx)))
            sparse_w=append(sparse_w,R[nzndx],1)
    # concatenate the i, j and w_ij into a single vector
    outlist=sparse_i
    outlist=append(outlist,sparse_j)
    outlist=append(outlist,sparse_w)
    # save the output file to a .NPY file
    save(outfile,outlist)
    print('finished ',infile,' len ',m)

### parcel_naming
"""
The following code was adapted from Satra Gosh's sad_figures.py script
located at https://github.com/satra/sad/blob/master/sad_figures.py
get cluster coordinates caculate the volume and the center of mass
of the brain regions in img
"""
def get_region_CoM(img, affine):
    coords = defaultdict(dict)
    # determine the unique regions
    labels = np.setdiff1d(np.unique(img.ravel()), [0])
    for label in labels:
        # calculate the volume of the region
        coords[label]["Vol"]=np.sum(img==label)
        # calculate the center of mass (CoM) of the region
        coords[label]["CoM"]=np.dot(affine,\
            np.hstack((np.mean(np.asarray(np.nonzero(img==label)),\
            axis = 1),1)))[:3].tolist()
    return (coords)

# convert coordinates from one image to another, coords'=inv(A2)*A1*coords
def most_common(lst):
    return max(set(lst), key=lst.count)

def image_downsample_voting(img, affine, down_img_template, down_img_affine):
    down_vals=defaultdict(list)
    old_coords=np.array(np.nonzero(img),dtype="int").T
    for i in range(np.shape(old_coords)[0]):
        new_coords=[str(int(c)) for c in np.round(\
            np.dot(np.linalg.inv(down_img_affine),\
            np.dot(affine,np.hstack((old_coords[i,:],1)))),decimals=0)[:3]]
        down_vals["_".join(new_coords)].append(img[tuple(old_coords[i,:])])
    new_img=np.zeros(np.shape(down_img_template),dtype="int")
    for k in down_vals.keys():
        idx=tuple([ int(n) for n in k.split("_")])
        new_img[idx]=most_common(down_vals[k])
    return (new_img)

def read_and_conform_atlas(atlas_file,atlas_label_file,\
        template_img,template_affine):
    atlas_labels=defaultdict()
    print("Reading in the atlas labels: %s"%(atlas_label_file))
    with open(atlas_label_file,"r") as f:
        for line in f:
            if '#' in line:
                continue
            line=line.rstrip('\n')
            vals=line.split(',')
            atlas_labels[int(vals[0])]=vals[1]
    atlas_labels[0]="None"
    print("Read in the atlas %s"%(atlas_file))
    # lets read in the Harvord Oxford Cortical map
    atlas_nii=nb.load(atlas_file)
    atlas_img=atlas_nii.get_data()
    print("Downsample the atlas")
    # resample the atlas to conform to parcels
    atlas_conform=image_downsample_voting(atlas_img, atlas_nii.get_affine(),\
                                   template_img, \
                                   template_affine);
    #print("Write out the downsampled atlas")
    #out_img=nb.Nifti1Image(atlas_conform,template_affine);
    #out_img.set_data_dtype("int16")
    #out_img.to_filename("atlas_conf.nii.gz")
    return (atlas_labels,atlas_conform)

def main():
    try:
        fsl_path=os.environ['FSLDIR']
    except KeyError:
        print("FSL_DIR is not set in the environment, is FSL installed?")
        sys.exit()
    # This is where the atlases are specified in the format
    # "ATLAS NAME":("path to label file","path to atlas.nii")
    atlas_cfg={\
    "Talairach Daemon":("talairach_labels.csv",\
      os.path.join(fsl_path,\
      "data/atlases/Talairach/Talairach-labels-2mm.nii.gz")),\
    "HarvardOxford Cortical":("harvardoxford_cortical_labels.csv",\
      os.path.join(fsl_path,\
      "data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz")),\
    "HarvardOxford Subcortical":("harvardoxford_subcortical_labels.csv",\
      os.path.join(fsl_path,\
      "data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz")),\
    "Jeulich Histological":("juelich_labels.csv",\
      os.path.join(fsl_path,\
      "data/atlases/Juelich/Juelich-maxprob-thr25-2mm.nii.gz")),\
    "MNI Structural":("mni_labels.csv",\
      os.path.join(fsl_path,\
      "data/atlases/MNI/MNI-maxprob-thr25-2mm.nii.gz"))}
    if len(sys.argv) < 4:
        print("number of arguements %d"%(len(sys.argv)))
        print("Usage %s <parcellation filename> <outname> <10,20,30,...>"%\
            (sys.argv[0]))
        sys.exit()
    parcel_filename=sys.argv[1]
    parcel_outname=sys.argv[2]
    parcel_vals=[int(n) for n in sys.argv[3].split(',')]
    print(parcel_vals)
    print("%s called with %s, %s, %s"%(sys.argv[0],parcel_filename,\
        parcel_outname,",".join([str(i) for i in parcel_vals])))
    print("Read in the parcellation results %s"%(parcel_filename))
    # lets read in the parcellation results that we want to label
    parcels_nii=nb.load(parcel_filename)
    parcels_img=parcels_nii.get_data()
    print(np.shape(parcels_img))
    print(len(parcel_vals))
    if len(parcel_vals) != np.shape(parcels_img)[3]:
        print("Length of parcel values (%d) != number of parcel images (%d)"%( \
            len(parcel_vals),np.shape(parcels_img)[3]))
        sys.exit()
    else:
        print("Length of parcel values (%d) == number of parcel images (%d)"%( \
            len(parcel_vals),np.shape(parcels_img)[3]))
    atlases=defaultdict()
    # read in the atlases
    for k in atlas_cfg.keys():
        atlases[k]=read_and_conform_atlas(atlas_cfg[k][1],atlas_cfg[k][0],\
            parcels_img[:,:,:,0],parcels_nii.get_affine())
    #for p in [0]:
    for p in range(np.shape(parcels_img)[3]):
        fid=open("%s_names_%d.csv"%(parcel_outname,parcel_vals[p]),"w")
        # print out the header
        fid.write("ROI number, volume, center of mass")
        for atlas in atlases:
            fid.write(",%s"%(atlas))
        fid.write("\n")
        p_c=get_region_CoM(parcels_img[:,:,:,p],parcels_nii.get_affine())
        for p_k in p_c.keys():
            fid.write("%d, %d, (%2.1f;%2.1f;%2.1f)"%(p_k,p_c[p_k]["Vol"],
                p_c[p_k]["CoM"][0],p_c[p_k]["CoM"][1],p_c[p_k]["CoM"][2]))
            for atlas in atlases.keys():
                fid.write(",")
                atlas_vals=atlases[atlas][1][np.nonzero(parcels_img[:,:,:,p]==p_k)]
                # calculate a histogram of the values
                atlas_hist=[(n,round(float(sum(atlas_vals==n))/float(len(atlas_vals)),2)) \
                    for n in np.unique(atlas_vals)]
                atlas_hist=sorted(atlas_hist,key=lambda f: -f[1])
                for h in atlas_hist:
                    if h[1] > 0.1:
                        fid.write("[\"%s\": %2.2f]"%(atlases[atlas][0][h[0]],h[1]))
            fid.write("\n")
        fid.close()

if __name__ == "__main__":
    main()

### python_ncut_lib
"""
(eigen_val, eigen_vec) = ncut( W, nbEigenValues ):

This function performs the first step of normalized cut spectral clustering.
The normalized LaPlacian is calculated on the similarity matrix W, and top
nbEigenValues eigenvectors are calculated. The number of eigenvectors
corresponds to the maximum number of classes (K) that will be produced by the
clustering algorithm.

   W:             symmetric #feature x #feature sparse matrix representing the
                  similarity between voxels, traditionally this matrix should
                  be positive semidefinite, but regularization is employed to
                  allow negative matrix entries (Yu 2001)
   nvEigenValues: number of eigenvectors that should be calculated, this
                  determines the maximum number of clusters (K) that can be
                  derived from the
   result
   eigen_val:     (output) eigenvalues from the eigen decomposition of the
                  LaPlacian of W
   eigen_vec:     (output) eigenvectors from the eign decomposition of the
                  LaPlacian of W
"""
def ncut( W, nbEigenValues ):
    # parameters
    offset=.5
    maxiterations=100
    eigsErrorTolerence=1e-6
    eps=2.2204e-16
    m=shape(W)[1]
    # make sure that W is symmetric, this is a computationally expensive
    # operation, only use for debugging
    # if (W-W.transpose()).sum() != 0:
    #    print "W should be symmetric!"
    #    exit(0)
    # degrees and regularization
    # S Yu Understanding Popout through Repulsion CVPR 2001
    # Allows negative values as well as improves invertability of d for small
    # numbers i bet that this is what improves the stability of the eigen
    d=abs(W).sum(0)
    dr=0.5*(d-W.sum(0))
    d=d+offset*2
    dr=dr+offset
    # calculation of the normalized LaPlacian
    W=W+spdiags(dr,[0],m,m,"csc")
    Dinvsqrt=spdiags((1.0/sqrt(d+eps)),[0],m,m,"csc")
    P=Dinvsqrt*(W*Dinvsqrt);
    # perform the eigen decomposition
    eigen_val,eigen_vec=eigsh(P,nbEigenValues,maxiter=maxiterations,\
        tol=eigsErrorTolerence,which='LA')
    # sort the eigen_vals so that the first
    # is the largest
    i=argsort(-eigen_val)
    eigen_val=eigen_val[i]
    eigen_vec=eigen_vec[:,i]
    # normalize the returned eigenvectors
    eigen_vec=Dinvsqrt*matrix(eigen_vec)
    norm_ones=norm(ones((m,1)))
    for i in range(0,shape(eigen_vec)[1]):
        eigen_vec[:,i]=(eigen_vec[:,i] / norm(eigen_vec[:,i]))*norm_ones
        if eigen_vec[0,i] != 0:
            eigen_vec[:,i] = -1 * eigen_vec[:,i] * sign( eigen_vec[0,i] )
    return(eigen_val, eigen_vec)

def discretisation( eigen_vec ):
    eps=2.2204e-16
    # normalize the eigenvectors
    [n,k]=shape(eigen_vec)
    vm=kron(ones((1,k)),sqrt(multiply(eigen_vec,eigen_vec).sum(1)))
    eigen_vec=divide(eigen_vec,vm)
    svd_restarts=0
    exitLoop=0
    ### if there is an exception we try to randomize and rerun SVD again
        ### do this 30 times
    while (svd_restarts < 30) and (exitLoop==0):
        # initialize algorithm with a random ordering of eigenvectors
        c=zeros((n,1))
        R=matrix(zeros((k,k)))
        R[:,0]=eigen_vec[int(rand(1)*(n-1)),:].transpose()
        for j in range(1,k):
            c=c+abs(eigen_vec*R[:,j-1])
            R[:,j]=eigen_vec[c.argmin(),:].transpose()
        lastObjectiveValue=0
        nbIterationsDiscretisation=0
        nbIterationsDiscretisationMax=20
        # iteratively rotate the discretised eigenvectors until they
        # are maximally similar to the input eignevectors, this
        # converges when the differences between the current solution
        # and the previous solution differs by less than eps or we
        # we have reached the maximum number of itarations
        while exitLoop == 0:
            nbIterationsDiscretisation = nbIterationsDiscretisation + 1
            # rotate the original eigen_vectors
            tDiscrete=eigen_vec*R
            # discretise the result by setting the max of each row=1 and
            # other values to 0
            j=reshape(asarray(tDiscrete.argmax(1)),n)
            eigenvec_discrete=csc_matrix((ones(len(j)),(range(0,n), \
                array(j))),shape=(n,k))
            # calculate a rotation to bring the discrete eigenvectors cluster to
            # the original eigenvectors
            tSVD=eigenvec_discrete.transpose()*eigen_vec
            # catch a SVD convergence error and restart
            try:
                U, S, Vh = svd(tSVD)
            except LinAlgError:
                # catch exception and go back to the beginning of the loop
                sys.stderr.write("SVD did not converge, randomizing and trying again")
                break
            # test for convergence
            NcutValue=2*(n-S.sum())
            if((abs(NcutValue-lastObjectiveValue) < eps ) or \
                      ( nbIterationsDiscretisation > \
                        nbIterationsDiscretisationMax )):
                exitLoop=1
            else:
                # otherwise calculate rotation and continue
                lastObjectiveValue=NcutValue
                R=matrix(Vh).transpose()*matrix(U).transpose()
    if exitLoop == 0:
        raise SVDError("SVD did not converge after 30 retries")
    else:
        return(eigenvec_discrete)
