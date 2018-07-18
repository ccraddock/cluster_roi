from unittest import TestCase


class TestConnectivity(TestCase):

    def test_cli_test_participant(self):
        import pyclusterroi as pc
        import os

        bids_dir = os.path.join(os.path.dirname(__file__), "test_data")

        # return_value = pc.pyclusterroi([bids_dir, bids_dir, 'test_config', '--skip_bids_validator', '-h'])
        # self.assertTrue(return_value == 0)

        return_value = pc.pyclusterroi(
            [bids_dir,
             bids_dir,
             'test_participant',
             '--skip_bids_validator',
             '--debug',
             '--mask', os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz',
             '--number_of_clusters', '10:100:10'])
        self.assertTrue(return_value == 0)

    def test_cli_participant(self):
        import pyclusterroi as pc
        import os

        bids_dir = os.path.join(os.path.dirname(__file__), "test_data")

        return_value = pc.pyclusterroi(
            [bids_dir,
             bids_dir,
             'participant',
             '--skip_bids_validator',
             '--mask', os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz',
             '--number_of_clusters', '10:100:10'])

        self.assertTrue(return_value == 0)

    def test_cli_group_2level(self):
        import pyclusterroi as pc
        import os

        bids_dir = os.path.join(os.path.dirname(__file__), "test_data")

        return_value = pc.pyclusterroi(
            [bids_dir,
             bids_dir,
             'group',
             '--group_name', 'allNFB',
             '--skip_bids_validator',
             '--mask', os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz',
             '--variant', 'tcorr90',
             '--number_of_clusters', '10:100:10'])
        self.assertTrue(return_value == 0)

    def test_cli_group_gmean(self):
        import pyclusterroi as pc
        import os

        bids_dir = os.path.join(os.path.dirname(__file__), "test_data")
        return_value = pc.pyclusterroi(
            [bids_dir,
             bids_dir,
             'group',
             '--group_name', 'allNFB',
             '--skip_bids_validator',
             '--group_level_method', 'gmean',
             '--mask', os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz',
             '--number_of_clusters', '10:100:10'])
        self.assertTrue(return_value == 0)

    def test_cli_participant_scorr_threads(self):
        import pyclusterroi as pc
        import os

        bids_dir = os.path.join(os.path.dirname(__file__), "test_data")

        return_value = pc.pyclusterroi(
            [bids_dir,
             bids_dir,
             'participant',
             '--skip_bids_validator',
             '--similarity_metric', 'scorr',
             '--thresh', '0.4',
             '--n_cpus', '2',
             '--mask', os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz',
             '--number_of_clusters', '10:100:10',
             '--participant_ndx', '1'])

        self.assertTrue(return_value == 0)

    # def test_make_ones_parcellation(self):
    #     """
    #     smoke test for ones connectivity
    #     :param self:
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import os
    #
    #     print('ones connectivity')
    #
    #     gm_mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #
    #     for participant in [1]:
    #
    #         participant_output_prefix = os.path.dirname(
    #             __file__) + "/test_data/derivatives/pyclusterroi/" + "sub-{0}/func/sub-{0}_task-rest_bold_".format(
    #             participant) + "variant-ones{k}_roi.nii.gz"
    #
    #         os.makedirs(os.path.dirname(participant_output_prefix), exist_ok=True)
    #
    #         out_files = pc.make_parcellation(range(10, 100, 10), 'ones', None, gm_mask_file,
    #                                          participant_output_prefix, thresh=0.5, num_threads=1)
    #
    #         self.assertTrue(out_files)
    #
    # def test_make_tcorr_parcellation(self):
    #     """
    #     smoke test for tcorr connectivity
    #     :param self:
    #     :return:
    #     """
    #
    #     import pyclusterroi as pc
    #     import os
    #
    #     gm_mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #
    #     print('tcorr connectivity')
    #
    #     for participant in [1]:
    #         participant_input_data = os.path.dirname(__file__) + \
    #                                  "/test_data/sub-{0}/func/sub-{0}_task-rest_bold.nii.gz".format(participant)
    #         participant_output_prefix = os.path.dirname(
    #             __file__) + "/test_data/derivatives/pyclusterroi/" + "sub-{0}/func/sub-{0}_task-rest_bold_".format(
    #             participant) + "variant-tcorr{k}_roi.nii.gz"
    #
    #         os.makedirs(os.path.dirname(participant_output_prefix), exist_ok=True)
    #
    #         out_files = pc.make_parcellation(range(10, 100, 10), 'tcorr', participant_input_data, gm_mask_file,
    #                                          participant_output_prefix, thresh=0.5, num_threads=1)
    #
    #         self.assertTrue(out_files)
    #
    # def test_make_tcorr_parcellation_mean(self):
    #     """
    #     smoke test for tcorr connectivity
    #     :param self:
    #     :return:
    #     """
    #
    #     import pyclusterroi as pc
    #     import os
    #
    #     print('tcorr connectivity')
    #
    #     gm_mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #
    #     participant_input_data = []
    #     for participant in [1, 2, 3]:
    #         participant_input_data.append(
    #             os.path.dirname(__file__) + "/test_data/sub-{0}/func/sub-{0}_task-rest_bold.nii.gz".format(participant))
    #
    #     participant_output_prefix = os.path.dirname(
    #         __file__) + "/test_data/derivatives/pyclusterroi/" + "group_task-rest_bold_variant-tcorr{k}_roi.nii.gz"
    #
    #     os.makedirs(os.path.dirname(participant_output_prefix), exist_ok=True)
    #
    #     out_files = pc.make_parcellation(range(10, 100, 10), 'tcorr', participant_input_data, gm_mask_file,
    #                                      participant_output_prefix, thresh=0.5, num_threads=1)
    #
    #     self.assertTrue(out_files)
    #
    # def test_make_tcorr_parcellation_2level(self):
    #     """
    #     smoke test for tcorr connectivity
    #     :param self:
    #     :return:
    #     """
    #
    #     import pyclusterroi as pc
    #     import os
    #     import multiprocessing as mp
    #
    #     class NoDaemonProcess(mp.Process):
    #         # make 'daemon' attribute always return False
    #         def _get_daemon(self):
    #             return False
    #
    #         def _set_daemon(self, value):
    #             pass
    #
    #         daemon = property(_get_daemon, _set_daemon)
    #
    #     # We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
    #     # because the latter is only a wrapper function, not a proper class.
    #     class MyPool(mp.pool.Pool):
    #         Process = NoDaemonProcess
    #
    #     print('tcorr connectivity')
    #
    #     gm_mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #     num_pt_at_once = 3
    #     num_threads = 2
    #
    #     with MyPool(num_pt_at_once) as pool:
    #         thread_results = []
    #         for participant in [1, 2, 3]:
    #             participant_input_data = os.path.dirname(
    #                 __file__) + "/test_data/sub-{0}/func/sub-{0}_task-rest_bold.nii.gz".format(participant)
    #
    #             participant_output_prefix = os.path.dirname(
    #                 __file__) + "/test_data/derivatives/pyclusterroi/" + "sub-{0}/func/sub-{0}_task-rest_bold_".format(
    #                 participant) + "variant-tcorr{k}_roi.nii.gz"
    #
    #             os.makedirs(os.path.dirname(participant_output_prefix), exist_ok=True)
    #
    #             thread_results.append(pool.apply_async(pc.make_parcellation, (100, 'tcorr', participant_input_data,
    #                                                                           gm_mask_file, participant_output_prefix,
    #                                                                           0.5, num_threads)))
    #
    #         participant_parcellation_data = []
    #         for thread_result in thread_results:
    #             participant_parcellation_data.append(thread_result.get()[0])
    #
    #     print(participant_parcellation_data)
    #
    #     group_output_prefix = os.path.dirname(
    #         __file__) + "/test_data/derivatives/pyclusterroi/group_task-rest_bold_variant-tcorr{k}_roi.nii.gz"
    #
    #     out_files = pc.make_parcellation(range(10, 100, 10), 'cluster', participant_parcellation_data, gm_mask_file,
    #                                      group_output_prefix, thresh=0.5, num_threads=1)
    #
    #     self.assertTrue(out_files)
    #
    # def test_make_scorr_parcellation(self):
    #     """
    #     smoke test for scorr connectivity
    #     :param self:
    #     :return:
    #     """
    #     import pyclusterroi as pc
    #     import os
    #
    #     gm_mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #
    #     print('scorr connectivity')
    #
    #     for participant in [1]:
    #         participant_input_data = os.path.dirname(__file__) + \
    #                                  "/test_data/sub-{0}/func/sub-{0}_task-rest_bold.nii.gz".format(participant)
    #         participant_output_prefix = os.path.dirname(
    #             __file__) + "/test_data/derivatives/pyclusterroi/" + "sub-{0}/func/sub-{0}_task-rest_bold_".format(
    #             participant) + "variant-scorr{k}_roi.nii.gz"
    #
    #         os.makedirs(os.path.dirname(participant_output_prefix), exist_ok=True)
    #
    #         ret_val = pc.make_parcellation(range(10, 100, 10), 'scorr', participant_input_data, gm_mask_file,
    #                                        participant_output_prefix, thresh=0.5, num_threads=8)
    #
    #         self.assertTrue(ret_val)
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
    #     import os
    #
    #     gm_mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #
    #     mask_array = nb.load(gm_mask_file).get_data()
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
    #     import os
    #
    #     gm_mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #
    #     mask_array = nb.load(gm_mask_file).get_data()
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
    #     import os
    #
    #     gm_mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #
    #     mask_array = nb.load(gm_mask_file).get_data()
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
    #     w_matrix = sp.csc_matrix((w, (i, j)), shape=(num_voxels, num_voxels))
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
    #     import os
    #
    #     print("testing make_local_connectivity_ones")
    #
    #     mask_file = os.path.dirname(__file__) + '/test_data/derivatives/gm_maskfile.nii.gz'
    #
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
    #     (w, i, j) = pc.make_local_connectivity("ones", None, mask_array, thresh=0.5, num_threads=1)
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
    #     import os
    #     import time
    #
    #     print("testing make_local_connectivity_tcorr")
    #
    #     mask_file = os.path.dirname(__file__) + "/test_data/derivatives/gm_maskfile.nii.gz"
    #     image_file = os.path.dirname(__file__) + "/test_data/sub-1/func/sub-1_task-rest_bold.nii.gz"
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
    #     start_time = time.time()
    #     (w, i, j) = pc.make_local_connectivity("tcorr", image_array, mask_array.reshape(mask_shape), thresh=0.5,
    #                                            num_threads=1)
    #     print("Calculated tcorr in {0} seconds".format(time.time() - start_time))
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
    #     import os
    #
    #     print("testing make_local_connectivity_scorr")
    #
    #     mask_file = os.path.dirname(__file__) + "/test_data/derivatives/gm_maskfile.nii.gz"
    #     image_file = os.path.dirname(__file__) + "/test_data/sub-1/func/sub-1_task-rest_bold.nii.gz"
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
    #     (w, i, j) = pc.make_local_connectivity('scorr', image_array, mask_array.reshape(mask_shape), thresh=0.5,
    #                                            num_threads=8)
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