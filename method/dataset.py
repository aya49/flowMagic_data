import os
import pandas as pd
import numpy as np

import torch
from torch.nn.functional import normalize
from torch.utils.data import Dataset

from transform import transform_dict

import random
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
        
        self.loadxy = opt.model !='setr'
        self.normx = False
        self.x_3D = False
        self.addclass = False
        self.ybig = False
        self.ysqueeze = False
        self.ymask = True
        self.addpos = False
        self.preload_data = opt.preload_data
        self.data_dir = opt.data_dir
        self.mode = opt.mode
        self.x_2D = opt.x_2D
        self.y_2D = opt.y_2D
        
        # seqlr = torch.zeros(opt.dim,opt.dim)
        # seqtb = torch.zeros(opt.dim,opt.dim)
        # for rci in range(0,opt.dim+1):
        #     seqlr[rci-1,:] = rci
        #     seqtb[:,rci-1] = rci
        # self.seqlr = 100*seqlr.float()/opt.dim
        # self.seqtb = 100*seqtb.float()/opt.dim
        
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
            self.x_contH = []
            self.x_contV = []
            self.y = []
            self.y_ = []
            goodi = []
            xyl = len(self.x_files[0])
            for i in range(xyl):
                try:
                    xil = []
                    for x2i in range(len(self.x_2D)):
                        loaded = pd.read_csv(self.x_files[x2i][i], header=None).values
                        if 'contourH' in self.x_2D[x2i]:
                            uH = np.unique(loaded)
                            uH.sort()
                        elif 'contourV' in self.x_2D[x2i]:
                            uV = np.unique(loaded)
                            uV.sort()
                        else:
                            xil.append(torch.tensor(loaded))
                    xil = torch.stack(xil)
                    
                    # yi = torch.tensor(pd.read_csv(self.y_files[i].replace(self.y_2D[0], '{}_'.format(self.y_2D[0])), header=None).values).unsqueeze(0)
                    yi = torch.tensor(pd.read_csv(self.y_files[i], header=None).values).unsqueeze(0)
                    
                    yi_ = torch.tensor(pd.read_csv(self.y_files[i].replace(self.y_2D[0], '{}_rough'.format(self.y_2D[0])), header=None).values).unsqueeze(0)
                    
                    self.x.append(xil)
                    self.x_contH.append(uH)
                    self.x_contV.append(uV)
                    self.y.append(yi)
                    self.y_.append(yi_)
                    goodi.append(i)
                    print(i, end=' ')
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
            self.y_ = torch.stack(self.y_)
        
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
    
    def factorize_labels(self):
        x_dirs_unique, x_dirs_factor = np.unique(np.array(self.x_dirs), return_inverse=True)
        self.x_dirs_factor = x_dirs_factor.tolist()
    
    def __getitem__(self, i):
        
        if self.preload_data:
            xi = self.x[i]
            if self.ymask:
                yi = self.y_[i]
            else:
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
            
            if self.ymask:
                yi = torch.tensor(pd.read_csv(self.y_files[i].replace(self.y_2D[0], '{}_rough'.format(self.y_2D[0])), header=None).values).unsqueeze(0)
            else:
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
        
        if self.transform != None:
            xi, yi = self.transform(xi, yi)
            
        if self.ysqueeze:
            yi = yi.squeeze()
        
        xi = xi.float()
        yi = yi.float()
        
        if len(self.x_contH)>0:
            thresH = self.x_contH[i]
            thresV = self.x_contV[i]
            cH = torch.ones(xi.shape[-2], xi.shape[-1])
            cV = torch.ones(xi.shape[-2], xi.shape[-1])
            if len(thresH)>1:
                for tH in range(1,len(thresH)):
                    cH[thresH[tH]-1:,:] = thresH[tH]
            if len(thresV)>1:
                for tV in range(1,len(thresV)):
                    cV[:,thresV[tV]-1:] = thresV[tV]
            
            chs = xi.shape[0]
            xit = [xi[ij] for ij in range(chs)]
            xit.extend([cH, cV])
            xi = torch.stack(xit)
        
        if self.addpos:
            chs = xi.shape[0]
            xit = [xi[ij] for ij in range(chs)]
            xit.extend([self.seqlr, self.seqtb])
            xi = torch.stack(xit)
        
        if self.normx:
            xi = normalize(xi)
        
        if self.x_3D:
            dimsize = xi.shape[1]
            xj = torch.zeros(1, dimsize, dimsize)
            xi = torch.cat((xi, xj))
        
        # experiment: add data/scat to data
        if self.addclass:
            xi[0][0][0] = self.x_dirs_factor[i]
        
        if self.ybig: # 3D y tensor
            yi = tensor2D3D(m=yi, C=self.n_class)
            
        if self.loadxy:
            return xi, yi
        
        return xi, yi, i, self.x_dirs[i], self.x_filenames[i]

# making this a separater function to split up Data2D into smaller chunks for saving
def split_Data2D(dataset, n):
    nsize = len(dataset) // n
    preload = dataset.preload_data
    
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

# get random subset of dataset
# making this a separater function to split up Data2D into smaller chunks for saving
def subset_Data2D(dataset, n):
    n_inds = random.sample(range(len(dataset)), n)
    preload = dataset.preload_data
    
    dataseti = copy.deepcopy(dataset)
    dataseti.x_dirs = [dataseti.x_dirs[i] for i in n_inds]
    dataseti.x_dirs_factor = [dataseti.x_dirs_factor[i] for i in n_inds]
    dataseti.x_filenames = [dataseti.x_filenames[i] for i in n_inds]
    if len(dataset.x_2D)>1:
        for j in range(len(dataset.x_2D)):
            dataseti.x_files[j] = [dataseti.x_files[j][i] for i in n_inds]
    else:
        dataseti.x_files = [dataseti.x_files[i] for i in n_inds]
    dataseti.y_files = [dataseti.y_files[i] for i in n_inds]
    
    if preload:
        dataseti.y_ = [dataseti.y_[i] for i in n_inds]
        dataseti.y = [dataseti.y[i] for i in n_inds]
        dataseti.x = [dataseti.x[i] for i in n_inds]
    
    return dataseti

def merge_Data2D(dataset, dataset_, preload=True):
    dataset.x_dirs.extend(dataset_.x_dirs)
    dataset.x_dirs_factor.extend(dataset_.x_dirs_factor)
    dataset.x_filenames.extend(dataset_.x_filenames)
    if len(dataset.x_2D) > 1:
        for j in range(len(dataset.x_2D)):
            dataset.x_files[j].extend(dataset_.x_files[j])
        else:
            dataset.x_files.extend(dataset_.x_files)
    dataset.y_files.extend(dataset_.y_files)
    
    if preload:
        dataset.y_ = torch.cat((dataset.y_, dataset_.y_), 0)
        dataset.y = torch.cat((dataset.y, dataset_.y), 0)
        dataset.x = torch.cat((dataset.x, dataset_.x), 0)
        dataset.x_contH.append(dataset_.x_contH)
        dataset.x_contV.append(dataset_.x_contV)
    return dataset


def tensor2D3D(m, C):
    m = m.squeeze()
    H, W = m.shape
    # C = int(torch.max(m))
    C_ = int(torch.max(m))
    r = torch.zeros(C, H, W)
    for i in range(C_+1):
        r[i][m==i] = 1
    
    return r

def tensor2D3D_(m, C):
    m = m.squeeze()
    N, H, W = m.shape
    # C = int(torch.max(m))
    C_ = int(torch.max(m))
    r = torch.zeros(N, C, H, W)
    for i in range(C_+1):
        for j in range(N):
            r[j][i][m[j]==i] = 1
    
    return r
