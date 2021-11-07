# set directory
import os
os.chdir("/home/aya43/flowMagic_data/src/method")
# os.chdir("/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src/method")

# basic modules
import sys
import time
import random

import copy # for deep copy
import pandas as pd
import numpy as np
import compress_pickle # pickle, but compresses
# import imp # update module with imp.reload()

# status modules
import trace
import inspect # inspect.getfullargspec(function) or signature # get function arguments
from GPUtil import showUtilization as gpu_usage # gpu_usage()

from matplotlib import pyplot as plt
# print to file
# orig_stdout = sys.stdout
# f = open('model_setr.txt', 'w')
# sys.stdout = f
# OUT
# sys.stdout = orig_stdout
# f.close()

# torch
import torch
from torch.utils.data import DataLoader
# from torchviz import make_dot # creates an image of the model
from torchsampler import ImbalancedDatasetSampler as ids # pip install https://github.com/ufoym/imbalanced-dataset-sampler/archive/master.zip

from opt import parse_options, update_opt
from util import prep_input, visualize
from transform import transform_dict
from dataset import Data2D, merge_Data2D, subset_Data2D
from models import create_model
from train_premeta import train
from util import load_checkpoint

## number of model parameters
# sum(p.numel() for p in model.parameters())
## If you want to calculate only the trainable parameters:
# sum(p.numel() for p in model.parameters() if p.requires_grad)

print("cuda available")
print(torch.cuda.is_available())


## OPTIONS ####################################################
opt = parse_options()
opt.preload_data = True # we pre-load everything so it's faster but takes up more memory
opt.num_workers = 32
opt.batch_size = 32 # if not enough gpu memory, reduce batch_size

# some convenient functions
def nomac(xx):
    return [x for x in xx if '__MACOSX' not in x]

def yegz(xx):
    return [x for x in xx if '.gz' in x]

flatx = lambda x: [i for row in x for i in row]


## DATA: datasets x 4 ########################################
dss = nomac( os.listdir(os.path.join(opt.data_folder, opt.x_2D[0])) )

x_dirs = nomac( flatx([[os.path.join(opt.data_folder, opt.x_2D[0], ds, sc) for 
                sc in os.listdir(os.path.join(opt.data_folder, opt.x_2D[0], ds))] for ds in dss]) )
x_dirs.sort()

# seqlr = torch.zeros(opt.dim,opt.dim)
# seqtb = torch.zeros(opt.dim,opt.dim)
# for rci in range(1, opt.dim+1):
#     seqlr[rci-1,:] = rci
#     seqtb[:,rci-1] = rci
# seqlr = 100*seqlr.float()/opt.dim
# seqtb = 100*seqtb.float()/opt.dim

# preload data
if opt.preload_data:
    ds_files = []
    for x_dir_mt in [x for x in x_dirs if 'pregnancy' in x]:
        xdmsplit = x_dir_mt.split('/')
        opt.data_scat = '/'.join(xdmsplit[-2:])
        ds_mt_r_path = os.path.join(opt.data_folder, 'dataloader_mt_r_{}.gz'.format(opt.data_scat.replace('/','_')))
        ds_files.append(ds_mt_r_path)
        # if not os.path.exists(ds_mt_r_path):
        print(x_dir_mt)
        x_files_mt = yegz(nomac( [os.path.join(x_dir_mt, f) for f in os.listdir(x_dir_mt)] ))
        dataset_mt_r = Data2D(opt, transform=transform_dict['A'], x_files=x_files_mt)
        compress_pickle.dump(dataset_mt_r, ds_mt_r_path, compression="lzma", set_default_extension=False) #gzip
