import time as time
from numpy import *
from scipy.sparse import csc_matrix
from python_ncut_lib import *


# group_mean_binfile_parcellate( infiles, outfile, K ):
#
# This function performs group level clustering of individual level clustering
# results. Each single subject clustering is converted into a coincidence
# matrix W, where w_ij = 1 if voxels i and j are in the same cluster, and zero
# otherwise. Coincidence matrices are averaged across subjects and then
# submitted to normalized cut clustering to obtain K clusters.
#    infiles:   list of .NPY or .bin file containing single subject clustering
#               results. Each of these contains a 3*N x 1 vector where the
#               first N values coorespond to the i indices, the next N
#               correspond to the j indices and the last N values correspond to
#               the weights w_ij. For more information on constructing these
#               files refer to make_local_connectivity_tcorr.py,
#               make_local_connectivity_scorr.py, or
#               make_local_connectivity_ones.py.
#    outfile:   a prefix for the output file, this name will be suffixed by
#               _K.npy where K corresponds to the clustering level
#    K:         list of numbers of clusters that will be generated. If this is
#               a single number then only that clustering will be generated. If
#               this is a list of numbers, then the normalized cut algorithm
#               will be run several times, once for each k in the list, and a
#               seperate output file will be generated for each clustering
#    n_voxels:  Number of voxels in the _mask_ used to generate the subject
#               specific connectivity matrices
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
