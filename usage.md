---
layout: page
title: "Usage Instructions and Tips"
description: "This is the usage page."
---
{% include JB/setup %}

The pyClusterROI_test.py script included in the distribution provides detailed
instructions and examples for using the library to parcellate fMRI data. This
script can be used in combination with a [**test
dataset**](https://www.nitrc.org/frs/downloadlink.php/3719) to validate your
installation. In addition to instructions for using the library functions,
there are several other details that bear consideration.

###Preparing data
Functional parcellation can be performed using either resting-state or
task-based fMRI data. The data should
be preprocessed to minimize the effects of confounding variation such as head
motion and physiological noise. Although the different preprocessing steps
applied to the data will have an impact on the clustering result, their impact
has yet to be systematically evaluated. The only mandatory preprocessing step
is that the data from different subjects must be spatially normalized to a
template space (e.g. MNI152 or Talaraich and Tournoux) so that brain regions
are aligned between them. The amount of spatial smoothing applied to the data
will also impact the clustering results, on one hand spatial smoothing improves
the correspondance between brain regions across individuals and improves
signal-to-noise ratio. On the other hand, too much spatial smoothing will bias
the clustering results. As with other preprocessing steps, the impact of
spatial smoothing on clustering results hasn't been evaluated, and the amount
of spatial smoothing applied should be considered carefully.

###Gray matter mask
A gray-matter mask is useful for reducing the computational complexity of
clustering, as well as, confining the clustering results to gray matter.
Creating the masks is  tricky, and the one that I used for the paper required
quite alot of trial and error. The process that I used was to segment each
individuals anatomical image and to average the GM probability maps across
subjects. After this I applied a threshold to the averaged probability maps in
order to create a binary GM mask. I chose a threshold by visual inspection so
that the mask covered all of the brain regions that I am interested in. You
will notice that parts of the putamen and globus pallidus may be excluded using
this procedure. You will have to try several different parameters for
segmentation and thresholding and possibly erosion/dilation in order to find a
happy medium. Another approach could be to use a brain atlas such as the
Harvard Oxford atlas to derive your grey matter mask. Of course, probably the
best way to go would be to hand segment a set of images ...

###Constructing connectivity matrices
Once the data has been preprocessed and a GM mask has been constructed, the
next step of functional parcellation is to construct individual level
similarity matrices from the data. This matrix is an Nvoxel x Nvoxel matrix in
which each entry corresponds to the similarity between voxels. Two different
methods are available for measuring similarity between voxels. Tcorr
corresponds to the Pearson's correlation between voxel time-series. Scorr is
the spatial correlation between whole-brain functional connectivity maps, each
of which are created by correlating the voxel time-series with every other
voxel time-series in the brain. The spatial constraint is imposed on the
clustering result by only calculating similarity between neigboring voxels. A
third method for generating connectivity matrices is included, the "ones
connectivity", sets the similarity between voxels to '1' if they are neighbors,
and '0' otherwise. In other words, only spatial information, and no functional
information, is used in the clustering.

The similarity metric chosen depends on the desired properties of the
clustering solution. Tcorr is preferred when the desire is to maximize the
temporal homogeniety of the clusters. Scorr should be used when the desire is
to maximaize the similarity between whole-brain functional connectivity maps
generated at each voxel. The spatial-only clustering is extremely fast and only
requires the gray-matter mask for construction. 

Instructions for generating individual-level connectivity matrices using each
of these methods are illustrated below (copined from pyClusterROI_test.py):

####Tcorr - temporal correlation

    {% highlight python %}
    from make_local_connectivity_tcorr import *

    # the name of the maskfile that we will be using
    maskname="gm_maskfile.nii.gz"

    # make a list of all of the input fMRI files that we 
    # will be using
    infiles = [  'subject1.nii.gz', 
                 'subject2.nii.gz', 
                 'subject3.nii.gz' ]

    # construct the connectivity matrices using tcorr 
    # and a r>0.5 threshold

    for idx, in_file in enumerate(infiles):

        # construct an output filename for this file
        outname='rm_tcorr_conn_'+str(idx)+'.npy'

        print 'tcorr connectivity',in_file
        # call the funtion to make connectivity
        make_local_connectivity_tcorr( in_file, 
            maskname, outname, 0.5 )
    {% endhighlight %}

####Scorr - spatial correlation

    {% highlight python %}
    from make_local_connectivity_scorr import *

    # the name of the maskfile that we will be using
    maskname="gm_maskfile.nii.gz"

    # make a list of all of the input fMRI files that we 
    # will be using
    infiles = [  'subject1.nii.gz',
                 'subject2.nii.gz',
                 'subject3.nii.gz' ]

    # construct the connectivity matrices using scorr 
    # and a r>0.5 threshold
    # This can take a _really_ long time
    for idx, in_file in enumerate(infiles):
    
        # construct an output filename for this file
        outname='rm_scorr_conn_'+str(idx)+'.npy'

        print 'scorr connectivity',in_file
        # call the funtion to make connectivity
        make_local_connectivity_scorr( in_file,
            maskname, outname, 0.5 )
    {% endhighlight %}

####Ones - spatial information only

    {% highlight python %}
    from make_local_connectivity_ones import *

    # the name of the maskfile that we will be using
    maskname="gm_maskfile.nii.gz"

    # the easiest is random clustering which doesn't 
    # require any functional data, just the mask
    print 'ones connectivity'
    make_local_connectivity_ones( maskname,
        'rm_ones_connectivity.npy')
    {% endhighlight %}


###Group-mean and two-level clustering
Clustering can begin once the individual level similarity matrices have been
constructed. pyClusterROI includes two different methods for combining
information across individuals to accomplish group-level clustering. In the
group-mean approach, the individual similarity matrices are averaged and the
resulting matrix is clustered using the normalized-cut algorithm. The two-level
approach begins with the normalized-cut clustering of each individuals
similarity matrix. Results are combined across individuals to for a
"coincidence matrix", that is clustered to obtain the group-level clustering
results. The "ones" clustering is accomplished by a single normalized-cut
clustering of the ones connectivity matrix. Below are examples for each of
these clustering methods (copied from pyClusterROI_test.py).

####Group-mean clustering

    {% highlight python %}
    from group_mean_binfile_parcellation import *
    NUM_CLUSTERS = [100,150,200]

    # group_mean clustering is pretty simple, input the 
    # connectivity files and run. We can perform multiple
    # clusterings with this function, so once again the
    # output filename is a prefix
    tcorr_conn_files=['rm_tcorr_conn_0.npy',
                      'rm_tcorr_conn_1.npy',
                      'rm_tcorr_conn_2.npy']

    print 'group-mean parcellate tcorr'
    group_mean_binfile_parcellate( tcorr_conn_files,\
        'rm_group_mean_tcorr_cluster', NUM_CLUSTERS,mask_voxels);
    {% endhighlight %}

####Two-level clustering

    {% highlight python %}
    from binfile_parcellation import *
    from group_binfile_parcellation import *

    NUM_CLUSTERS = [100,150,200]

    # for tcorr
    for idx, in_file in enumerate(infiles):

        # construct filenames
        infile='rm_tcorr_conn_'+str(idx)+'.npy'
        outfile='rm_tcorr_indiv_cluster_'+str(idx)

        print 'tcorr parcellate',in_file
        binfile_parcellate(infile, outfile, NUM_CLUSTERS)

    # the 2-level clustering has to be performed once for 
    # each desired clustering level, and requires individual 
    # level clusterings as inputs
    for k in NUM_CLUSTERS:
        ind_clust_files=[]
        for i in range(0,len(infiles)):
            ind_clust_files.append('rm_tcorr_indiv_cluster_'+\
	        str(i)+'_'+str(k)+'.npy')

        print '2-level parcellate tcorr',k
        group_binfile_parcellate(ind_clust_files,\
            'rm_group_tcorr_cluster_'+str(k)+'.npy',k,mask_voxels)
    {% endhighlight %}


####Ones clustering

    {% highlight python %}
    from binfile_parcellation import *

    NUM_CLUSTERS = [100,150,200]

    # For random custering, this is all we need to do, there is no 
    # need for group level clustering, remember that the output 
    # filename is a prefix
    binfile_parcellate('rm_ones_connectivity.npy',\
        'rm_ones_cluster',NUM_CLUSTERS)
    {% endhighlight %}
    
####Finalizing the cluster atlas
The outputs of the clustering must be transformed in order to make them more
easy to use. The output .npy files are 'index files' and must be remapped to
voxel space. Additionally, the regions are not consecutively numbered. Utility
functions are provided to handle these transformations.

    {% highlight python %}
    from make_image_from_bin_renum import *
    # write out for group mean clustering
    for k in NUM_CLUSTERS:
        binfile='rm_group_mean_tcorr_cluster_'+str(k)+'.npy'
        imgfile='rm_group_mean_tcorr_cluster_'+str(k)+'.nii.gz'
        make_image_from_bin_renum(imgfile,binfile,maskname)
    {% endhighlight %}


####Labelling clusters
Several individuals have asked for labels for the clusters. A python script is
now provided that calculates the center of mass for each cluster and maps this
coordinate to brain region labels. This mapping is performed for several
different brain atlases, including Harvard-Oxford, Juelich, Talaraich and
Tournoux, and the MNI atlases. This script is called from the command line
using the following syntax:

    {% highlight python %}
    python parcel_naming.py tcorr05_2level_all.nii.gz \
        tcorr05_2level '50,100,150,200'
    {% endhighlight %}
