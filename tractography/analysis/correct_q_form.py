import nibabel as nib
import numpy as np


path = "subj1_tensor.nii.gz"

tensor = nib.load(path)
tensor.set_qform(tensor.get_qform(), 1)

nib.save(tensor, path)
