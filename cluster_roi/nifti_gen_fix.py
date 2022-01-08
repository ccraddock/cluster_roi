import sys
import numpy as np
import nibabel as nib

"""
 Temporary band-aid script to generate NIfTI images from NumPY parcellation data.

 Written by Dan Lurie (danjlurie@gmail.com).

 Usage:

 python nifti_gen_fix.py /path/to/parcel_data.npy /path/to/mask.nii.gz /path/to/new_atlas_image.nii.gz

 """

 # Read in command line arguments.
 parcel_path, mask_path, out_path = sys.argv[1:4]

 # Load parcellation data.
 parcel_data = np.load(parcel_path)

 # Load the mask file.
 mask_file = nib.load(mask_path)

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

 # Write the NIfTI file.
 atlas_image.to_filename(out_path)
