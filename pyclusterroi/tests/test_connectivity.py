from unittest import TestCase

#import nibabel as nb




class TestConnectivity(TestCase):

    #smoke test for bin connectivity
    def test_make_bin_connectivity(self):
        import pyclusterroi as pc
        import nibabel as nb
        import numpy as np
        print 'ones connectivity'

        maskarray=nb.load("../../gm_maskfile.nii.gz").get_data()

        mvals = pc.make_local_connectivity_ones(maskarray)

        np.save("../../ones_connectivity.npy",mvals)

        self.assertTrue(isinstance(mvals, list))

    #smoke test for tcorr connectivity
    def test_make_tcorr_connectivity(self):
        import pyclusterroi as pc
        import nibabel as nb
        import numpy as np
        print 'tcorr connectivity'

        maskarray=nb.load("../../gm_maskfile.nii.gz").get_data()
        imdat=nb.load("../../subject1.nii.gz").get_data()

        mvals = pc.make_local_connectivity_tcorr(imdat,maskarray,0.6)
        print np.shape(mvals)

        self.assertTrue(isinstance(mvals, list))

    #smoke test for scorr connectivity
    def test_make_scorr_connectivity(self):
        import pyclusterroi as pc
        import nibabel as nb
        import numpy as np
        print 'scorr connectivity'

        maskarray=nb.load("../../gm_maskfile.nii.gz").get_data()
        imdat=nb.load("../../subject1.nii.gz").get_data()

        mvals = pc.make_local_connectivity_scorr(imdat,maskarray,0.6)
        print np.shape(mvals)

        self.assertTrue(isinstance(mvals, list))


    def test_get_neighbors_3d(self):
        import pyclusterroi as pc
        import random
        import numpy as np
        import nibabel as nb

        mask_array=nb.load("../../gm_maskfile.nii.gz").get_data()
        mask_size=np.shape(mask_array)

        ix=np.nonzero(mask_array.flatten())[0]

        test_ix = ix[random.sample(range(0,len(ix)),1)][0]

        print mask_size

        d_ndx = pc.get_neighbors(test_ix, mask_size)
        print d_ndx
        assert (len(d_ndx) == 26)


    def test_get_neighbors_3d(self):
        import pyclusterroi as pc
        import random
        import numpy as np
        import nibabel as nb

        mask_array = nb.load("../../gm_maskfile.nii.gz").get_data()
        mask_size = np.shape(mask_array)

        ix = np.nonzero(mask_array.flatten())[0]

        test_ix = ix[random.sample(range(0, len(ix)), 1)][0]
        d_ndx = pc.get_neighbors(test_ix, mask_size)

        assert (len(d_ndx) == 26)

    def test_get_neighbors_2d(self):
        import pyclusterroi as pc
        import random
        import numpy as np
        import nibabel as nb

        mask_array = nb.load("../../gm_maskfile.nii.gz").get_data()
        mask_array = mask_array[24,:,:]
        mask_size = np.shape(mask_array)

        ix = np.nonzero(mask_array.flatten())[0]

        test_ix = ix[random.sample(range(0, len(ix)), 1)][0]

        d_ndx = pc.get_neighbors(test_ix, mask_size)
        assert (len(d_ndx) == 8)

    def test_get_neighbors_1d(self):
        import pyclusterroi as pc
        import random
        import numpy as np
        import nibabel as nb

        mask_array = nb.load("../../gm_maskfile.nii.gz").get_data()
        mask_array = mask_array[24,28,:]
        mask_size = np.shape(mask_array)

        ix = np.nonzero(mask_array.flatten())[0]

        test_ix = ix[random.sample(range(0, len(ix)), 1)][0]

        d_ndx = pc.get_neighbors(test_ix, mask_size)
        assert (len(d_ndx) == 2)
