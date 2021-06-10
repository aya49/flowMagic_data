## DataLoader class loads 2D samples
## author: alice yue aya43@sfu.ca
## date created: 2021-03-21
## date modified: 2021-03-21

import os
import pandas as pd
import numpy as np
import torch as pt
from torch.utils.data import Dataset, DataLoader

# dataset loader class: torch.utils.data.Dataset
# overwrite class methods:
# - __len__ so that len(dataset) returns the size of the dataset. READ CSV HERE
# - __getitem__ to support indexing such that dataset[i] can be used to get i th sample. READ IMAGES HERE
# our data set is a dict {'image': image, 'landmarks': landmarks} + optional argument transform for any preprocessing
class DataLoader(Dataset):
  
  def __init__(self, data_dir, x_dens=None, x_2Ddensity="x_2Ddensity", x_2Dscatter="x_2Dscatter", y_2D="y_2D"):
    """
    args:
      data_dir: 400 * 400 dataset directory (file format: .csv.gz) containing the following folders:
        x_2Ddensity (smoothed 2D Gaussian kernel density)
        x_2Dscatter (1/0 whether or not there is a cell on pixel)
        y_2D (numeric label matrix; 0 means other)
        (not used; constant removed in derivative) y_2Dncells (numeric count matrix with number of cells on each pixel)
      x_dens: None; otherwise it is a list of x_2Ddensity directories
      * folder structure: data_dir > x_2Ddensity/x_2Dscatter/y_2D > dataset > scatterplot > sample.csv.gz
    """
    
    # lambda to flatten lists without importing functions...
    flatx = lambda x: [i for row in x for i in row]
    
    # list x_dens data set directories and files
    if x_dens is None:
      x_dens = flatx([[os.path.join(data_dir, x_2Ddensity, ds, sc)
                       for sc in os.listdir(os.path.join(data_dir, x_2Ddensity, ds))]
                      for ds in os.listdir(os.path.join(data_dir, x_2Ddensity))])
    
    x_dens_files = flatx([[os.path.join(x_den, f) for f in os.listdir(x_den)] for x_den in x_dens])
    x_scat_files = [x_dens_file.replace(x_2Ddensity, x_2Dscatter) for x_dens_file in x_dens_files]
    y_2D_files = [x_dens_file.replace(x_2Ddensity, y_2D) for x_dens_file in x_dens_files]
    
    # self.variables
    self.data_dir = data_dir          # data set root directory
    self.dscats = [os.path.dirname(x_dens_file).replace(os.path.join(data_dir, x_2Ddensity), "")[1:] for x_dens_file in x_dens_files]
    self.x_dens_files = x_dens_files  # data set files: x_2Ddensity
    self.x_scat_files = x_scat_files  # data set files: x_2Dscatter
    self.y_2D_files = y_2D_files      # data set files: y_2D
  
  def __len__(self):
    return len(self.xdens)
  
  def __getitem__(self, i):
    xdi = pt.tensor(pd.read_csv(self.x_dens_files[i], header=None).values)
    xsi = pt.tensor(pd.read_csv(self.x_scat_files[i], header=None).values)
    xi = pt.cat([xdi.unsqueeze_(-1), xsi.unsqueeze_(-1)], dim=2).squeeze_()
    yi = pt.tensor(pd.read_csv(self.y_2D_files[i], header=None).values)
    
    return xi, yi, self.dscats[i]