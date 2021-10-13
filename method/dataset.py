import os
import pandas as pd
import numpy as np

import torch
from torch.utils.data import Dataset

from transform import transform_dict

import copy

# dataset loader class: torch.utils.data.Dataset
# overwrite class methods:
# - __len__ so that len(dataset) returns the size of the dataset. READ CSV HERE
# - __getitem__ to support indexing such that dataset[i] can be used to get i th sample. READ IMAGES HERE
# our data set is a dict {'image': image, 'landmarks': landmarks} + optional argument transform for any preprocessing
class Data2D(Dataset):
  
    def __init__(self, opt, transform=None, x_dirs=None, x_files=None):
        """
        args:
        data_dir: 200 * 200 dataset directory (file format: .csv.gz) containing the following folders:
            x_2D (list of input data directories; 1st directory corresponds to x_dirs)
                x_2Ddenscat (smoothed 2D Gaussian kernel density; 0 on pixels where there are no cells)
                x_2Dcontour (0-1 1=line density contours)
            y_2D (numeric label matrix; 0 means other)
            (not used; constant removed in derivative) y_2Dncells (numeric count matrix with number of cells on each pixel)
        x_dirs: None; or list of x_2Ddenscat directories (pretrain), list of train files (metatrain), list of test files (metatest)
        mode: "pretrain", "distill", "meta" (meta="train" or "test")
        * folder structure: data_dir > x_2Ddenscat/x_2Dcontour/y_2D > dataset > scatterplot > sample.csv.gz
        """
        
        # lambda to flatten lists without importing functions...
        flatx = lambda x: [i for row in x for i in row]
        
        # list x_dirs data set directories and files
        if x_dirs is not None:
            if type(x_dirs) != list:
                x_dirs = [x_dirs]
        
        # if 'meta' in opt.mode: # x_dirs is just one directory
        #     x_files_ = [os.path.join(opt.x_2D[0], opt.data_scat, fn) for 
        #                     fn in os.listdir(opt.shot_dir)]
        #     x_files_ = [x for x in x_files_ if '__MACOSX' not in x]
        #             
        #     if opt.mode == 'metatest':
        #         # x_files = list of test denscat samples of 1 scatterplot
        #         x_files = [os.path.join(opt.x_2D[0], opt.data_scat, fn) for 
        #                     fn in os.listdir(os.path.join(opt.x_2D[0], opt.data_scat))]
        #         x_files = [x for x in x_files if x not in x_files_ and '__MACOSX' not in x]
        #     else:
        #         # x_files = list of reference denscat samples of 1 scatterpot
        #         x_files = x_files_
        #     
        # else:
        if x_files is None:
            if x_dirs is None:
                x_dirs = flatx([[os.path.join(opt.data_dir, opt.x_2D[0], ds, sc) for 
                            sc in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0], ds))] for 
                            ds in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0]))])
                x_dirs = [x for x in x_dirs if '__MACOSX' not in x]
            x_files = flatx([flatx([[os.path.join(x_den, f) for f in os.listdir(x_den)] for x_den in x_dirs])])
            x_files = [x for x in x_files if '__MACOSX' not in x]
        elif x_dirs is None:
            x_dirs = [os.path.dirname(x_file) for x_file in x_files]
        
        # factorize data/scatterplot
        x_dirs_unique, x_dirs_factor = np.unique(np.array(x_dirs), return_inverse=True)
        x_dirs_factor = x_dirs_factor.tolist()

        y_files = [x_file.replace(opt.x_2D[0], opt.y_2D[0]) for x_file in x_files]

        if len(opt.x_2D) > 1:
            x_files = [x_files]
            for x2i in range(1, len(opt.x_2D)):
                x_files.append([x_file.replace(opt.x_2D[0], opt.x_2D[x2i]) for x_file in x_files[0]])

        self.preload_data = opt.preload_data
        self.data_dir = opt.data_dir          # data set root directory
        self.mode = opt.mode
        self.x_2D = opt.x_2D
        self.y_2D = opt.y_2D

        # if self.mode == 'metatest':
        #     self.ycell = list([])
        #     self.ydiscrete_files = [x_file.replace(opt.x_2D[0], opt.y_2D[1]) for x_file in x_files[0]]
        #    self.yvector_files = [x_file.replace(opt.x_2D[0], opt.y_2D[2]) for x_file in x_files[0]]

        self.transform = transform
        self.n_class = opt.n_class

        # these could be  @property; def x_dirs ...
        self.x_dirs = x_dirs              # scatterplot ("class")
        self.x_dirs_factor = x_dirs_factor
        self.x_filenames = [os.path.basename(x) for x in x_files[0]]
        self.x_files = x_files            # list of data set files: x_2D[0]
        self.y_files = y_files            # data set files: y_2D

        self.data_dir = opt.data_dir

        if opt.preload_data: # or len(x_files) < 200 # *** change
            self.x = []
            self.y = []
            goodi = []
            xyl = len(self.x_files[0])
            for i in range(xyl):
                print(i)
                try:
                    xil = []
                    for x2i in range(len(self.x_2D)):
                        xil.append(torch.tensor(pd.read_csv(self.x_files[x2i][i].replace(self.x_2D[x2i], self.x_2D[0]), header=None).values))
                    xil = torch.stack(xil)
                    self.x.append(xil)

                    yi = torch.tensor(pd.read_csv(self.y_files[i], header=None).values).unsqueeze(0)
                    self.y.append(yi)
                    goodi.append(i)
                except:
                    print("error")
                    pass
            
            x_dirs = self.x_dirs
            x_dirs_factor = self.x_dirs_factor
            x_filenames = self.x_filenames
            x_files = self.x_files
            y_files = self.y_files

            if isinstance(goodi, int):
                goodi = [goodi]
            if (len(range(xyl))>len(goodi)):
                self.x_dirs = [x_dirs[i] for i in goodi]
                self.x_dirs_factor = [x_dirs_factor[i] for i in goodi]
                self.x_filenames = [x_filenames[i] for i in goodi]
                if (len(self.x_2D)>1):
                    for j in range(len(x_files)):
                        self.x_files[j] = [x_files[j][i] for i in goodi]
                else:
                    self.x_files = [x_files[i] for i in goodi]
                self.y_files = [y_files[i] for i in goodi]

            # if round(i/xyl, 2) == prog:
            #     print('{prog} ')
            #     prog += .05

            self.x = torch.stack(self.x)
            self.y = torch.stack(self.y)
            
        
        # if 'meta' in self.mode:
        #     self.ydiscrete = list([])
        #     for i in range(len(self.ydiscrete_files)):
        #         yi0 = torch.tensor(pd.read_csv(self.ydiscrete_files[i], header=None).values)
        #         yi = yi0.squeeze()
        #         self.ydiscrete.append(yi)
        # 
        #     self.yvector = list([])
        #     for i in range(len(self.yvector_files)):
        #         yi0 = torch.tensor(pd.read_csv(self.yvector_files[i], header=None).values)
        #         yi = yi0.squeeze()
        #         self.yvector.append(yi)


    def __len__(self):
        return len(self.x_files[0])
    
    def get_labels(self):
        return self.x_dirs_factor
    
    def transform(self, transform):
        self.transform = transform
    
    def __getitem__(self, i):

        if self.preload_data:
            xi = self.x[i]
            yi = self.y[i]
            # if 'meta' in self.mode:
            #     ydi = self.ydiscrete[i]
            #     yvi = self.yvector[i]

        else:
            xi = []
            for x2i in range(len(self.x_2D)):
                xi.append(torch.tensor(pd.read_csv(self.x_files[x2i][i].replace(self.x_2D[x2i], self.x_2D[0]), header=None).values))
            xi = torch.stack(xi)
            # xi = torch.cat(xil, dim=0).squeeze_()
            
            yi = torch.tensor(pd.read_csv(self.y_files[i], header=None).values).unsqueeze(0)
            # yi = torch.zeros(self.n_class, yi0.shape[0], yi0.shape[1])
            # for yc in range(torch.max(i).int()):
            #     yi[yc, yi0 == (yc + 1)] = 1

            # if 'meta' in self.mode:
            #     ydi0 = torch.tensor(pd.read_csv(self.ydiscrete_files[i], header=None).values)
            #     ydi = ydi0.squeeze()
            #
            #     yvi0 = torch.tensor(pd.read_csv(self.yvector_files[i], header=None).values)
            #     yvi = yvi0.squeeze()
        
        # if 'meta' in self.mode:
        #     return xi, yi, i, self.x_dirs[i], ydi, yvi

        # experiment: add data/scat to data
        xi[0][0][0] = self.x_dirs_factor[i]

        if self.transform != None:
            xi, yi = self.transform(xi, yi)
        # yi = yi.squeeze()
        
        return xi, yi, i, self.x_dirs[i], self.x_filenames[i]

# making this a separater function to split up Data2D into smaller chunks for saving
def split_Data2D(dataset, n, preload=True):
    nsize = len(dataset) // n
    
    datasets = []
    for i in range(n):
        if i == (n-1):
            datasets.append(dataset)
        else:
            dataseti = copy.deepcopy(dataset)
            
            dataseti.x_dirs = dataseti.x_dirs[0:nsize]
            dataset.x_dirs = dataset.x_dirs[nsize:]
            dataseti.x_dirs_factor = dataseti.x_dirs_factor[0:nsize]
            dataset.x_dirs_factor = dataset.x_dirs_factor[nsize:]
            dataseti.x_filenames = dataseti.x_filenames[0:nsize]
            dataset.x_filenames = dataset.x_filenames[nsize:]
            if len(dataset.x_2D)>1:
                for j in range(len(dataset.x_2D)):
                    dataseti.x_files[j] = dataseti.x_files[j][0:nsize]
                    dataset.x_files[j] = dataset.x_files[j][nsize:]
            else:
                dataseti.x_files = dataseti.x_files[0:nsize]
                dataset.x_files = dataset.x_files[nsize:]
            dataseti.y_files = dataseti.y_files[0:nsize]
            dataset.y_files = dataset.y_files[nsize:]
            
            if preload:
                dataseti.y = dataseti.y[0:nsize]
                dataset.y = dataset.y[nsize:]
                dataseti.x = dataseti.x[0:nsize]
                dataset.x = dataset.x[nsize:]
            datasets.append(dataseti)
    return datasets

def merge_Data2D(datasets, preload=True):
    dataset = copy.deepcopy(datasets[0])
    for i in range(1, len(datasets)):
        dataset = datasets[i]
        dataset.x_dirs.append(datasets[i].x_dirs)
        dataset.x_dirs_factor.append(datasets[i].x_dirs_factor)
        dataset.x_filenames.append(datasets[i].x_filenames)
        if len(dataset.x_2D) > 1:
            for j in range(len(dataset.x_2D)):
                dataset.x_files[j].append(datasets[i].x_files[j])
            else:
                dataset.x_files.append(datasets[i].x_files)
        dataset.y_files.append(datasets[i].y_files)
        
        if preload:
            dataset.y = torch.cat((dataset.y, datasets[i].y), 0)
            dataset.x = torch.cat((dataset.x, datasets[i].x), 0)
    return dataset


