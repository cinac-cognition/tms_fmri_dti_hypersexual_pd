
import glob
import os
from pathlib import Path
import pydicom


alldcm = glob.glob('./data/original_dicom/*/*/*.dcm')
for jj in range(0,len(alldcm)):
    ds = pydicom.dcmread(alldcm[jj])
    if jj is 0:
        studyUID = ds.StudyInstanceUID
    ds.StudyInstanceUID= studyUID
    ds.save_as(alldcm[jj])

for i in range(0, len(alldcm)):
    ds = pydicom.dcmread(alldcm[i])
    fn = os.path.basename(alldcm[i])
    path = Path('/path/to/symlinks/' + ds.StudyInstanceUID)
    path.mkdir(parents=True, exist_ok=True)
    print(alldcm[i], ds.StudyInstanceUID)
    os.symlink(alldcm[i], '/path/to/symlinks/' + ds.StudyInstanceUID + '/' + fn)
