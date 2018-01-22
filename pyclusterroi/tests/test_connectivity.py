from unittest import TestCase


class TestConnectivity(TestCase):


    def test_make_ones_parcellation(self):
        """
        smoke test for ones connectivity
        :param self:
        :return:
        """
        import pyclusterroi as pc
        print 'ones connectivity'

        for pt in [1,2,3]:
            retVal = pc.make_individual_parcellation(range(10,100,10),'ones',
                                                     '../../test_in/sub-{pt}/func/sub-{pt}_task-rest_bold.nii.gz'.format(pt),
                                                     '../../test_in/derivatives/gm_maskfile.nii.gz'.format(pt),
                                                     '../../test_in/derivatives/sub-{pt}/func/spatncut/sub-{pt}_task-rest_bold_var-ones\{nroi\}_roi.nii.gz'.format(pt),
                                                     thresh=0.5)

        self.assertTrue(retVal)

    # def test_make_tcorr_parcellation(self):
    #     """
    #     smoke test for tcorr connectivity
    #     :param self:
    #     :return:
    #     """
    #
    #     import pyclusterroi as pc
    #
    #     print 'tcorr connectivity'
    #
    #     ret_val = pc.make_individual_parcellation(range(10, 100, 10), 'tcorr', '../../subject1.nii.gz',
    #                                               '../../gm_maskfile.nii.gz', 'tcorr.nii.gz', thresh=0.5)
    #
    #     self.assertTrue(ret_val)
    #
    # def test_make_scorr_parcellation(self):
    #     """
    #     smoke test for scorr connectivity
    #     :param self:
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #
    #     print 'scorr connectivity'
    #
    #     ret_val = pc.make_individual_parcellation(range(10, 100, 10), 'scorr', '../../subject1.nii.gz',
    #                                               '../../gm_maskfile.nii.gz', 'scorr.nii.gz', thresh=0.5)
    #
    #     self.assertTrue(ret_val)
    #
    # def test_get_neighbors_3d(self):
    #     """
    #     smoke for get neighbors functionality
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import random
    #     import numpy as np
    #     import nibabel as nb
    #
    #     mask_array = nb.load("../../gm_maskfile.nii.gz").get_data()
    #     mask_size = np.shape(mask_array)
    #
    #     ix = np.nonzero(mask_array.flatten())[0]
    #
    #     test_ix = ix[random.sample(range(0, len(ix)), 1)][0]
    #
    #     d_ndx = pc.get_neighbors(test_ix, mask_size)
    #
    #     assert (len(d_ndx) == 26)
    #
    # def test_get_neighbors_2d(self):
    #     """
    #     smoke for get 2d neighbors functionality
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import random
    #     import numpy as np
    #     import nibabel as nb
    #
    #     mask_array = nb.load("../../gm_maskfile.nii.gz").get_data()
    #     mask_array = mask_array[24, :, :]
    #     mask_size = np.shape(mask_array)
    #
    #     ix = np.nonzero(mask_array.flatten())[0]
    #
    #     test_ix = ix[random.sample(range(0, len(ix)), 1)][0]
    #
    #     d_ndx = pc.get_neighbors(test_ix, mask_size)
    #     assert (len(d_ndx) == 8)
    #
    # def test_get_neighbors_1d(self):
    #     """
    #     smoke test for test_get_neighbors_1d functionality
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import random
    #     import numpy as np
    #     import nibabel as nb
    #
    #     mask_array = nb.load("../../gm_maskfile.nii.gz").get_data()
    #     mask_array = mask_array[24, 28, :]
    #     mask_size = np.shape(mask_array)
    #
    #     ix = np.nonzero(mask_array.flatten())[0]
    #
    #     test_ix = ix[random.sample(range(0, len(ix)), 1)][0]
    #
    #     d_ndx = pc.get_neighbors(test_ix, mask_size)
    #     assert (len(d_ndx) == 2)
    #
    # def test_make_local_connectivity_clusters(self):
    #     """
    #     smoke test for test_get_neighbors_1d functionality
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import numpy as np
    #     import scipy.sparse as sp
    #
    #     print("testing make_local_connectivity_clusters")
    #
    #     num_clusters = 100
    #     num_voxels = 18343
    #
    #     im_array = np.random.choice(num_clusters, num_voxels, replace=True)
    #
    #     w, i, j = pc.make_local_connectivity_clusters(im_array)
    #
    #     w_matrix = sp.csc_matrix((w,(i,j)), shape=(num_voxels, num_voxels))
    #
    #     self.assertTrue(np.mean(w) == 1.0)
    #     self.assertTrue(np.min(i) >= 0 and np.max(i) < num_voxels)
    #     self.assertTrue(np.min(j) >= 0 and np.max(j) < num_voxels)
    #     self.assertTrue((w_matrix - w_matrix.transpose()).nnz == 0)
    #
    # def test_make_local_connectivity_ones(self):
    #     """
    #     smoke test for test_get_neighbors_1d functionality
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import numpy as np
    #     import nibabel as nb
    #     import scipy.sparse as sp
    #
    #     print("testing make_local_connectivity_ones")
    #
    #     mask_file = "../../gm_maskfile.nii.gz"
    #     # load mask
    #     nim = nb.load(mask_file)
    #     mask_array = ((nim.get_data()) > 0).astype('uint32')
    #
    #     # make the mask do joint duty as a lookup table
    #     num_mask_vx = 0
    #     for i in range(0, mask_array.shape[0]):
    #         for j in range(0, mask_array.shape[1]):
    #             for k in range(0, mask_array.shape[2]):
    #                 if mask_array[i, j, k] == 1:
    #                     num_mask_vx += 1
    #                     mask_array[i, j, k] = num_mask_vx
    #
    #     w, i, j = pc.make_local_connectivity_ones(mask_array)
    #
    #     w_matrix = sp.csc_matrix((w, (i, j)), shape=(num_mask_vx, num_mask_vx))
    #
    #     self.assertTrue(np.mean(w) == 1.0)
    #     self.assertTrue(np.min(i) >= 0 and np.max(i) < num_mask_vx)
    #     self.assertTrue(np.min(j) >= 0 and np.max(j) < num_mask_vx)
    #     self.assertTrue((w_matrix - w_matrix.transpose()).nnz == 0)
    #
    # def test_make_local_connectivity_tcorr(self):
    #     """
    #     smoke test for test_get_neighbors_1d functionality
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import numpy as np
    #     import nibabel as nb
    #     import scipy.sparse as sp
    #
    #     print("testing make_local_connectivity_tcorr")
    #
    #     mask_file = "../../gm_maskfile.nii.gz"
    #     image_file = "../../subject1.nii.gz"
    #
    #     # load mask
    #     nim = nb.load(mask_file)
    #     mask_array = ((nim.get_data()) > 0).astype('uint32')
    #     mask_array = mask_array[12:25, 14:42, 12:34]
    #
    #     mask_shape = mask_array.shape
    #     mask_array = mask_array.flatten()
    #
    #     # make the mask do joint duty as a lookup table
    #     num_mask_voxels = 0
    #     for mask_index, mask_value in enumerate(mask_array):
    #         if mask_value == 1:
    #             num_mask_voxels += 1
    #             mask_array[mask_index] = num_mask_voxels
    #
    #     print("number of in-mask voxels: {0}".format(num_mask_voxels))
    #
    #     # load and conform the fMRI data
    #     image_array = nb.load(image_file).get_data()
    #     image_array = image_array[12:25, 14:42, 12:34]
    #
    #     image_shape = image_array.shape
    #     num_time_courses = image_shape[-1]
    #
    #     # reshape fmri data to a num_voxels x num_time_points array
    #     image_array = np.reshape(image_array, (np.prod(image_shape[:-1]), num_time_courses))
    #
    #     # reduce fmri data to just in-brain voxels
    #     image_array = image_array[mask_array != 0, :]
    #
    #     (w, i, j) = pc.make_local_connectivity_tcorr(image_array, mask_array.reshape(mask_shape), thresh=0.5)
    #
    #     w_matrix = sp.csc_matrix((w, (i, j)), shape=(num_mask_voxels, num_mask_voxels))
    #
    #     self.assertTrue(np.min(w) >= 0.5 and np.max(w) <= 1.0)
    #     self.assertTrue(np.min(i) >= 0 and np.max(i) < num_mask_voxels)
    #     self.assertTrue(np.min(j) >= 0 and np.max(j) < num_mask_voxels)
    #     self.assertTrue((w_matrix - w_matrix.transpose()).nnz == 0)
    #
    # def test_make_local_connectivity_scorr(self):
    #     """
    #     smoke test for make_local_connectivity_scorr functionality
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import numpy as np
    #     import nibabel as nb
    #     import scipy.sparse as sp
    #
    #     print("testing make_local_connectivity_scorr")
    #
    #     mask_file = "../../gm_maskfile.nii.gz"
    #     image_file = "../../subject1.nii.gz"
    #
    #     # load mask
    #     nim = nb.load(mask_file)
    #     mask_array = ((nim.get_data()) > 0).astype('uint32')
    #     mask_array = mask_array[12:25, 14:42, 12:34]
    #
    #     mask_shape = mask_array.shape
    #     mask_array = mask_array.flatten()
    #
    #     # make the mask do joint duty as a lookup table
    #     num_mask_voxels = 0
    #     for mask_index, mask_value in enumerate(mask_array):
    #         if mask_value == 1:
    #             num_mask_voxels += 1
    #             mask_array[mask_index] = num_mask_voxels
    #
    #     print("number of in-mask voxels: {0}".format(num_mask_voxels))
    #
    #     # load and conform the fMRI data
    #     image_array = nb.load(image_file).get_data()
    #     image_array = image_array[12:25, 14:42, 12:34]
    #
    #     image_shape = image_array.shape
    #     num_time_courses = image_shape[-1]
    #
    #     # reshape fmri data to a num_voxels x num_time_points array
    #     image_array = np.reshape(image_array, (np.prod(image_shape[:-1]), num_time_courses))
    #
    #     # reduce fmri data to just in-brain voxels
    #     image_array = image_array[mask_array != 0, :]
    #
    #     (w, i, j) = pc.make_local_connectivity_scorr(image_array, mask_array.reshape(mask_shape), thresh=0.5)
    #
    #     print("size of w {0}, i {1}, j {2}".format(len(w), len(i), len(j)))
    #
    #     w_matrix = sp.csc_matrix((w, (i, j)), shape=(num_mask_voxels, num_mask_voxels))
    #
    #     self.assertTrue(np.min(w) >= 0.5 and np.max(w) <= 1.0)
    #     self.assertTrue(np.min(i) >= 0 and np.max(i) < num_mask_voxels)
    #     self.assertTrue(np.min(j) >= 0 and np.max(j) < num_mask_voxels)
    #
    #     # for some reason there appears to be more quantization error with scorr than tcorr
    #     # lets call anything withing 1e-6 close enough
    #     self.assertTrue(((w_matrix-w_matrix.transpose()).data < 1e-6).all())