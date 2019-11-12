import os
import numpy as np
import nrrd

str_path = 'experiment_562061175' # Your ID number
str_path_right = str_path +'_right'
loadpath = r'I:\python_code\mouse_connectivity_models-master\paper\connectivity' # Your folder path
os.makedirs(os.path.join(loadpath, str_path_right))
fpath = os.path.join(loadpath, str_path ,'projection_density_100.nrrd')
PROJ, metaPROJ = nrrd.read(fpath)
PROJ = PROJ[..., ::-1]
nrrd.write(os.path.join(loadpath, str_path_right, 'projection_density_100.nrrd'),PROJ)

fpath = os.path.join(loadpath, str_path ,'data_mask_100.nrrd')
PROJ, metaPROJ = nrrd.read(fpath)
PROJ = PROJ[..., ::-1]
nrrd.write(os.path.join(loadpath, str_path_right, 'data_mask_100.nrrd'),PROJ)

fpath = os.path.join(loadpath, str_path ,'injection_density_100.nrrd')
PROJ, metaPROJ = nrrd.read(fpath)
PROJ = PROJ[..., ::-1]
nrrd.write(os.path.join(loadpath, str_path_right, 'injection_density_100.nrrd'),PROJ)

fpath = os.path.join(loadpath, str_path ,'injection_fraction_100.nrrd')
PROJ, metaPROJ = nrrd.read(fpath)
PROJ = PROJ[..., ::-1]
nrrd.write(os.path.join(loadpath, str_path_right, 'injection_fraction_100.nrrd'),PROJ)

# Change folder name
IDname = os.path.join(loadpath, str_path)
str_path_left = str_path +'_left'
newname = os.path.join(loadpath, str_path_left)
os.rename(IDname, newname)
oldname = os.path.join(loadpath, str_path_right)
os.rename(oldname, IDname)