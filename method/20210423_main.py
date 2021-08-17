# module load LIB/CUDA/10.2

# !pip install segmentation_models_pytorch
# !pip install torchviz
# !MMCV_WITH_OPS=1 pip install mmcv-full==1.2.2 -f https://download.openmmlab.com/mmcv/dist/cu110/torch1.7.0/index.html #full version not for windows

import sys
import os
import trace

import pandas as pd
import numpy as np

import torch
from torch.utils.data import Dataset
import torchvision.transforms as transforms
import torchvision.transforms.functional as FT
import torch.nn as nn
import torch.nn.functional as FN
from torch.utils.data import DataLoader

from torchviz import make_dot
import segmentation_models_pytorch as smp

import time
import random

import mmcv
import mmseg

from mmcv.utils import Config
from mmseg.models import build_segmentor

from opt import parse_options
from transform import transform_dict
from dataset import Data2D
from datasampler import ImbalancedDatasetSampler as ids
from torch.utils.data import DataLoader

from model import create_model

from train import train

print("cuda available")
print(torch.cuda.is_available())

# options
opt = parse_options()

# tensorboard logger
logger = tb_logger.Logger(logdir=opt.tb_folder, flush_secs=2)


## DATA: datasets x 4 ###########################################
dss = os.listdir(os.path.join(opt.data_dir, opt.x_2D[0]))

# PARAMETER
dti = 0

ds_tr = dss[-dti]
ds_mt = dss[dti]

# list of denscat folders
flatx = lambda x: [i for row in x for i in row]
x_dirs_tr = flatx([[os.path.join(opt.data_dir, opt.x_2D[0], ds, sc) for 
             sc in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0], ds))] for 
             ds in ds_mt])
random.shuffle(x_dirs_tr)
val_len = round(len(x_dirs_tr)/10)
x_dirs_tr_tr = x_dirs_tr[val_len:]
x_dirs_tr_val = x_dirs_tr[:val_len]

x_dirs_mts = [os.path.join(opt.data_dir, opt.x_2D[0], ds_mt, sc) for 
             sc in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0], ds_mt))]



# ## DATA: scatterplot style x   #################################
# dss = os.listdir(os.path.join(opt.data_dir, opt.x_2D[0]))

# # list of denscat folders
# flatx = lambda x: [i for row in x for i in row]
# x_dirs = flatx([[os.path.join(opt.data_dir, opt.x_2D[0], ds, sc) for 
#          sc in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0], ds))] for 
#          ds in dss])

# # PARAMETER (only lowest scoring one in each style will act as meta)
# sti = [3,5,6,12,20]
# smi = 0
# # dti = [3,5,6,12,20]
# # dmi = 0

# x_dirs_tr = x_dirs[sti]
# x_dirs_mt = x_dirs[smi]


## PRE TRAIN #################################################
if opt.mode == 'pretrain':
    # create dataloader
    dataset_tr = Data2D(opt, transform=transform_dict['A'], x_dirs=x_dirs_tr_tr)
    dataloader_tr = DataLoader(sampler=ids(dataset_tr), batch_size=opt.batch_size,
                    shuffle=True, drop_last=True, num_workers=opt.num_workers)

    dataset_val = Data2D(opt, transform=transform_dict['A'], x_dirs=x_dirs_tr_val)
    dataloader_val = DataLoader(dataset_val)

    # initialize model
    model = create_model(opt)

    # train
    train(opt, model, dataloader_tr, dataloader_val) # opt.preload_model = True


# ## DISTILL: work in progress ##################################
# if opt.mode == 'distill':
#     # load model
#     model   = create_model(opt.model, opt.n_class, opt.dim)
#     model_t = create_model(opt.model, opt.n_class, opt.dim)
#     ckpt = torch.load(opt.model_path, map_location="cuda:0" if torch.cuda.is_available() else "cpu")
#     model_t.load_state_dict(torch.load(ckpt)['model'])

#     train(opt, model, dataloader_tr, model_t)
    


for n_shot in [1, 2, 3, 4, 5, 10, 15, 20]:
    for x_dirs_mts_files in x_dirs_mts:
        ## META TRAIN ################################################
        # load model
        model = create_model(opt.model, opt.n_class, opt.dim)
        ckpt = torch.load(opt.model_path, map_location="cuda:0" if torch.cuda.is_available() else "cpu")
        model.load_state_dict(torch.load(ckpt)['model'])

        # get n-shot samples
        mt_dir = os.path.join(x_dirs_mts_files[0].split('/')[0:-1])
        pam_dir = mt_dir.replace("data/2D/x_2Ddenscat","data/2D/x_2Ddenscat_euclidean_rankkmed")
        mt_files = os.path.join(mt_dir, os.listdir(pam_dir + "/" + str(n_shot)))
        mv_files = [mvf for mvf in x_dirs_mts_files if mvf not in mt_files]
        
        # create dataloader
        dataset_tr = Data2D(opt, transform=transform_dict['B'], x_dirs=mt_files)
        dataloader_tr = DataLoader(sampler=ids(dataset_tr), batch_size=opt.batch_size,
                        shuffle=True, drop_last=True, num_workers=opt.num_workers)

        dataset_val = Data2D(opt, x_dirs=mt_files)
        dataloader_val = DataLoader(dataset_val)

        train(opt, model, dataloader_tr, dataloader_val)

        ## META TEST #################################################
        start = time.time()
        val_acc, val_std = meta_test(model, meta_valloader)
        print('val_acc: {:.4f}, val_std: {:.4f}, time: {:.1f}'.format(val_acc, val_std, time.time() - start))

        start = time.time()
        val_acc_feat, val_std_feat = meta_test(model, meta_valloader, use_logit=False)
        val_time = time.time() - start
        print('val_acc_feat: {:.4f}, val_std: {:.4f}, time: {:.1f}'.format(val_acc_feat, val_std_feat, time.time() - start))

        start = time.time()
        test_acc, test_std = meta_test(model, meta_testloader)
        print('test_acc: {:.4f}, test_std: {:.4f}, time: {:.1f}'.format(test_acc, test_std, time.time() - start))

        start = time.time()
        test_acc_feat, test_std_feat = meta_test(model, meta_testloader, use_logit=False)
        print('test_acc_feat: {:.4f}, test_std: {:.4f}, time: {:.1f}'.format(test_acc_feat, test_std_feat, time.time() - start))


