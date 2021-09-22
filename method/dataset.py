import os
import pandas as pd

import torch
from torch.utils.data import Dataset

import pickle


# dataset loader class: torch.utils.data.Dataset
# overwrite class methods:
# - __len__ so that len(dataset) returns the size of the dataset. READ CSV HERE
# - __getitem__ to support indexing such that dataset[i] can be used to get i th sample. READ IMAGES HERE
# our data set is a dict {'image': image, 'landmarks': landmarks} + optional argument transform for any preprocessing
class Data2D(Dataset):
  
    def __init__(self, opt, transform=None, x_dirs=None):
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
        
        if 'meta' in opt.mode: # x_dirs is just one directory
            x_files_ = [os.path.join(opt.x_2D[0], opt.data_scat, fn) for 
                            fn in os.listdir(opt.shot_dir)]
                    
            if opt.mode == 'metatest':
                # x_files = list of test denscat samples of 1 scatterplot
                x_files = [os.path.join(opt.x_2D[0], opt.data_scat, fn) for 
                            fn in os.listdir(os.path.join(opt.x_2D[0], opt.data_scat))]
                x_files = [x for x in x_files if x not in x_files_]
            else:
                # x_files = list of reference denscat samples of 1 scatterpot
                x_files = x_files_
            
        else:
            if x_dirs is None:
                x_dirs = flatx([[os.path.join(opt.data_dir, opt.x_2D[0], ds, sc) for 
                            sc in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0], ds))] for 
                            ds in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0]))])
            x_files = [flatx([[os.path.join(x_den, f) for 
                            f in os.listdir(x_den)] for 
                            x_den in x_dirs])]
        
        
        y_files = [x_file.replace(opt.x_2D[0], opt.y_2D[0]) for x_file in x_files[0]]

        if (len(opt.x_2D) > 0):
            for x2i in range(1, len(opt.x_2D)):
                x_files.append([x_file.replace(opt.x_2D[0], opt.x_2D[x2i]) for x_file in x_files[0]])
        
        x_dirs = [os.path.dirname(x_file) for x_file in x_files]

        self.preload_data = opt.preload_data
        self.data_dir = opt.data_dir          # data set root directory
        self.mode = opt.mode
        self.meta = opt.meta
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
        self.x_files = x_files            # list of data set files: x_2D[0]
        self.y_files = y_files            # data set files: y_2D

        self.data_dir = opt.data_dir

        # self.preload_data = opt.preload_data # or len(x_files) < 200 # *** change
        if self.preload_data:
            self.x = list([])
            for i in range(len(x_files)):
                xil = list([])
                for x2i in range(len(self.x_2D)):
                    xil.append(torch.tensor(pd.read_csv(self.x_files[i].replace(self.x_2D[0], self.x_2D[x2i]), header=None).values).unsqueeze_(0))
                xi = torch.cat(xil, dim=0).squeeze_()
                self.x.append(xi)
            
            self.y = list([])
            for i in range(len(y_files)):
                yi0 = torch.tensor(pd.read_csv(self.y_files[i], header=None).values)
                yi = yi0.unsqueeze(0)
                self.y.append(yi)
            
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
        return len(self.x_files)
    
    def pickle(self):
        f = file(os.path.join(self.data_dir, 'data_temp'), 'wb')
        pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def unpickle(self):
        with file(os.path.join(self.data_dir, 'data_temp'), 'rb') as f:
            return pickle.load(f)
    
    def __getitem__(self, i):

        if self.preload_data:
            xi = self.x[i]
            yi = self.y[i]
            # if 'meta' in self.mode:
            #     ydi = self.ydiscrete[i]
            #     yvi = self.yvector[i]

        else:
            xil = list([])
            for x2i in range(len(self.x_2D)):
                xil.append(torch.tensor(pd.read_csv(self.x_files[i].replace(self.x_2D[0], self.x_2D[x2i]), header=None).values).unsqueeze_(0))
            xi = torch.cat(xil, dim=0).squeeze_()
            
            yi0 = torch.tensor(pd.read_csv(self.y_files[i], header=None).values)
            yi = yi0.unsqueeze(0)
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

        if self.transform != None:
            xi, yi = self.transform(xi, yi, xi.shape[-1])
        yi = yi.squeeze()
        
        return xi, yi, i, self.x_dirs[i]


