import nibabel as nb
import numpy as np
from collections import defaultdict


# The following code was adapted from Satra Gosh's sad_figures.py script
# located at https://github.com/satra/sad/blob/master/sad_figures.py

# get cluster coordinates caculate the volume and the center of mass
# of the brain regions in img

def get_coords(img, affine):
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

