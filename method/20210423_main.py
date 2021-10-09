# this training was done on 1 gpu and 32 workers

# for image classification: https://github.com/WangYueFt/rfs

import os
os.chdir("/home/aya43/flowMagic_data/src/method")
# os.chdir("/mnt/FCS_local3/backup/Brinkman group/current/Alice/flowMagic_data/src/method")
import inspect # inspect.getfullargspec(function) or signature # get function arguments

import sys
import trace

import pandas as pd
import numpy as np

import random

import torch
from torch.utils.data import Dataset
import torchvision.transforms as transforms
import torchvision.transforms.functional as FT
import torch.nn as nn
import torch.nn.functional as FN
from torch.utils.data import DataLoader

# from torchviz import make_dot

import time
import random

import mmcv
import mmseg

import compress_pickle

from mmcv.utils import Config
# from mmseg.models import build_segmentor

from opt import parse_options
from transform import transform_dict
from dataset import Data2D
# from datasampler import ImbalancedDatasetSampler as ids
from torchsampler import ImbalancedDatasetSampler as ids
# pip install https://github.com/ufoym/imbalanced-dataset-sampler/archive/master.zip
from torch.utils.data import DataLoader

# import imp # update module with imp.reload()

## number of parameters
# sum(p.numel() for p in model.parameters())
## If you want to calculate only the trainable parameters:
# sum(p.numel() for p in model.parameters() if p.requires_grad)

import tensorboard_logger as tb_logger

from GPUtil import showUtilization as gpu_usage # gpu_usage()

from models import create_model

from train import train, validate

print("cuda available")
print(torch.cuda.is_available())

# print to file
# orig_stdout = sys.stdout
# f = open('model_setr.txt', 'w')
# sys.stdout = f
# OUT
# sys.stdout = orig_stdout
# f.close()

# options
opt = parse_options()

# tensorboard logger
# logger = tb_logger.Logger(logdir=opt.tb_folder, flush_secs=2)

# I used a mac, so I need to get rid of '__MACOSX'
def nomac(xx):
    return [x for x in xx if '__MACOSX' not in x]


## DATA: datasets x 4 ########################################
dss = nomac( os.listdir(os.path.join(opt.data_dir, opt.x_2D[0])) )

## PRE-TRAIN #################################################
# choose the data set 0-3 we use as the metatest data set
## for dti in range(4):
dti = 3 ##
ds_tr = [x for i, x in enumerate(dss) if i!=dti]
ds_mt = dss[dti]

# train/metatrain data sets denscats folder paths
flatx = lambda x: [i for row in x for i in row]
x_dirs_tr = nomac( flatx([[os.path.join(opt.data_dir, opt.x_2D[0], ds, sc) for 
            sc in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0], ds))] for ds in ds_tr]) )
x_dirs_mt = nomac( [os.path.join(opt.data_dir, opt.x_2D[0], ds_mt, sc) for 
            sc in os.listdir(os.path.join(opt.data_dir, opt.x_2D[0], ds_mt))] )

## if opt.mode == 'pretrain':
# split pre-train data set into train (95%) and validation (5%)
x_files_tr = nomac( flatx([flatx([[os.path.join(x_den, f) for f in os.listdir(x_den)] for x_den in x_dirs_tr])]) )

x_files_tr_v_ind = random.sample(range(0,len(x_files_tr)), int(len(x_files_tr)/20))
x_files_tr_v = [x_files_tr[x] for x in range(0,len(x_files_tr)) if x in x_files_tr_v_ind]
x_files_tr_t = [x_files_tr[x] for x in range(0,len(x_files_tr)) if x not in x_files_tr_v_ind]

# set some parameters --- if not enough gpu memory, reduce batch_size
opt.preload_data = True # we pre-load everything so it's faster but takes up more memory
opt.num_workers = 32
opt.batch_size = 32
opt.cuda = 'cuda:0'

# create datasets
ds_tr_t_path = os.path.join(opt.data_dir, 'dataset_tr_t_{}.gz'.format(ds_mt))
if os.path.exists(ds_tr_t_path):
    dataset_tr_t = compress_pickle.load(ds_tr_t_path, compression="lzma", set_default_extension=False) #gzip
else:
    dataset_tr_t = Data2D(opt, transform=transform_dict['A'], x_files=x_files_tr_t)
    compress_pickle.dump(dataset_tr_t, ds_tr_t_path, compression="lzma", set_default_extension=False) #gzip

ds_tr_v_path = os.path.join(opt.data_dir, 'dataset_tr_v_{}.gz'.format(ds_mt))
if os.path.exists(ds_tr_v_path):
    dataset_tr_v = compress_pickle.load(ds_tr_v_path, compression="lzma", set_default_extension=False) #gzip
else:
    dataset_tr_v = Data2D(opt, transform=transform_dict['A'], x_files=x_files_tr_v)
    compress_pickle.dump(dataset_tr_v, ds_tr_v_path, compression="lzma", set_default_extension=False) #gzip

# create dataloaders
dataloader_tr_t = DataLoader(dataset=dataset_tr_t, sampler=ids(dataset_tr_t), 
                             batch_size=opt.batch_size, drop_last=True, #shuffle=True, 
                             num_workers=opt.num_workers)
dataloader_tr_v = DataLoader(dataset=dataset_tr_v,
                             batch_size=opt.batch_size, drop_last=False, shuffle=False,
                             num_workers=opt.num_workers)

# initialize model
model = create_model(opt).cuda()
# sum(p.numel() for p in model.parameters())

optimizer = torch.optim.Adam(model.parameters(), lr=opt.learning_rate, weight_decay=0.0005)

# train and validate
opt.epochs = 50
opt.save_freq = 10
acc, loss, model = train(opt=opt, model=model, train_loader=dataloader_tr_t, val_loader=dataloader_tr_v, optimizer=optimizer) # pt.preload_model = True


# ## DISTILL: work in progress ##################################
# if opt.mode == 'distill':
#     # load model
#     model   = create_model(opt.model, opt.n_class, opt.dim)
#     model_t = create_model(opt.model, opt.n_class, opt.dim)
#     ckpt = torch.load(opt.model_path, map_location="cuda:0" if torch.cuda.is_available() else "cpu")
#     model_t.load_state_dict(torch.load(ckpt)['model'])

#     train(opt, model, dataloader_tr, model_t)
        









## try just training with 10 samples ####################################
n_shot = 10
mt_files = []
mv_files = []
for x_dir_mt in x_dirs_mt:
    x_dirs_mts_files = [os.path.join(x_dir_mt, f) for f in os.listdir(x_dir_mt)]
    x_dirs_mts_files = [x for x in x_dirs_mts_files if '__MACOSX' not in x]
    # get n-shot samples
    mt_files.append(random.sample(x_dirs_mts_files, n_shot))
    mv_files.append([mvf for mvf in x_dirs_mts_files if mvf not in mt_files])

mt_files = flatx(mt_files)
mv_files = flatx(mv_files)

opt.save_dir = opt.save_dir + '_'
os.makedirs(opt.model_dir, exist_ok=True)
opt.model_folder = opt.model_folder + '_'
os.makedirs(opt.model_folder, exist_ok=True) 
opt.tb_dir = opt.tb_dir + '_'

dataset_tr_t = Data2D(opt, transform=transform_dict['A'], x_files=mt_files)
dataset_tr_v = Data2D(opt, transform=transform_dict['A'], x_files=mv_files)

opt.num_workers = 32
opt.batch_size = 16
opt.preload_data = True
opt.cuda = 'cuda:0'

dataloader_tr_t = DataLoader(dataset=dataset_tr_t, 
                    sampler=ids(dataset_tr_t), 
                    batch_size=opt.batch_size,# shuffle=True, 
                    drop_last=True, num_workers=opt.num_workers)
dataloader_tr_v = DataLoader(dataset=dataset_tr_v,
                    batch_size=opt.batch_size // 2, shuffle=False, drop_last=False,
                    num_workers=opt.num_workers // 2)

# initialize model
model = create_model(opt).cuda(device=opt.cuda)
# sum(p.numel() for p in model.parameters())

optimizer = torch.optim.Adam(model.parameters(), lr=opt.learning_rate, weight_decay=0.0005)

# train and validate
opt.epochs = 5000
opt.save_freq = 50
acc, loss, model = train(opt=opt, model=model, train_loader=dataloader_tr_t, val_loader=dataloader_tr_v, optimizer=optimizer) # pt.preload_model = True

# get results
dataset_tr_v.transform = transform_dict['B']
dataloader_tr_test = DataLoader(dataset=dataset_tr_v,
                    batch_size=len(dataset_tr_v), shuffle=False, drop_last=False,
                    num_workers=1)
for idx, (inp, target, i, xdir, xfn) in enumerate(dataloader_tr_test):
    break

set_cuda = torch.cuda.is_available()

inp = inp.float()
inp = inp.cuda(device=opt.cuda) if set_cuda else inp
target = target.cuda(device=opt.cuda) if set_cuda else target


model.eval()

start_i = 0
xdir_ = xdir[0]
(H, W, C) = (256, 256, len(opt.x_2D))

acc = []

for j in range(len(inp)-1):
    if xdir[j+1] != xdir_:
        end_i = j
        if j == len(inp)-1:
            inp_ = inp[:][start_i:(j+1)]
            target_ = target[:][start_i:(j+1)]
            img_metas_ = [{
                'img_shape': (H, W, C),
                'ori_shape': (H, W, C),
                'pad_shape': (H, W, C),
                'filename': xfn__,
                'scale_factor': 1.0,
                'flip': False,
            } for xfn__ in xfn[start_i:]]
        else:
            inp_ = inp[:][start_i:]
            target_ = target[:][start_i:]
            img_metas_ = [{
                'img_shape': (H, W, C),
                'ori_shape': (H, W, C),
                'pad_shape': (H, W, C),
                'filename': xfn__,
                'scale_factor': 1.0,
                'flip': False,
            } for xfn__ in xfn[start_i:]]
            xdir_ = xdir[j+1]
        
        scores = model.forward(inp_, img_metas_, gt_semantic_seg=target_, return_loss=True)
        acc.append([xdir_, float(scores['decode.acc_seg'])])
        # res = model.inference(inp, img_metas, rescale=False)

        if j != len(inp)-1:
            xdir_ = xdir[j+1]
        start_i = j+1






## if opt.mode == 'meta':
opt.mode = 'meta'
## for n_shots in [1, 2, 3, 4, 5, 10, 15, 20]:
n_shots = 10
opt.n_shots = n_shots
## for x_dir_mt in x_dirs_mt:
x_dir_mt = x_dirs_mt[0]
xdmsplit = x_dir_mt.split('/')
opt.data_scat = '/'.join(xdmsplit[-2:])

x_files_mt = nomac( [os.path.join(x_dir_mt, f) for f in os.listdir(x_dir_mt)] )

# get n-shot samples
opt.model_name_meta = '{}_METAdatascat:{}_METAshots:{}'.format(opt.model_name, opt.data_scat, opt.n_shots) # data_scat e.g. 'pregnancy/07_FoxP3CD25_CD4Tcell'
opt.shot_dir = os.path.join(opt.shot_dir, opt.data_scat, str(opt.n_shots) + '.csv.gz')
x_files_mt_t = pd.read_csv(opt.shot_dir)
# file handling HERE!!!
x_files_mt_t =  random.sample(x_files_mt, opt.n_shots) ## TEMP!!!!

# get test samples
x_files_mt_r = list(set(x_files_mt) - set(x_files_mt_t))


## META-TRAIN #################################################
dataset_mt_r = Data2D(opt, transform=transform_dict['B'], x_files=x_files_mt_r)

# set some parameters --- if not enough gpu memory, reduce batch_size
opt.num_workers = 32
opt.batch_size = 32
opt.preload_data = True # we pre-load everything so it's faster but takes up more memory
opt.cuda = 'cuda:0'

# create datasets
dataset_mt_t = Data2D(opt, transform=transform_dict['A'], x_files=x_files_mt_t)
dataset_mt_v = dataset_mt_t
dataset_mt_v.transform = transform_dict['B']

# create dataloaders
dataloader_mt_t = DataLoader(dataset=dataset_mt_t, sampler=ids(dataset_mt_t), 
                             batch_size=opt.batch_size, drop_last=True, # shuffle=True, 
                             num_workers=opt.num_workers)
dataloader_mt_v = DataLoader(dataset=dataset_mt_v, sampler=ids(dataset_mt_v), 
                             batch_size=len(dataset_mt_v), drop_last=False, shuffle=False, 
                             num_workers=opt.num_workers)

# load model
model = create_model(opt).cuda()
ckpt = torch.load(os.path.join(opt.model_folder, '{}_last.pth'.format(opt.model)))
model.load_state_dict(ckpt['model'])

optimizer = torch.optim.Adam(model.parameters(), lr=opt.learning_rate, weight_decay=0.0005)

# train and validate
opt.epochs = 50
opt.save_freq = 50
acc, loss, model = train(opt=opt, model=model, train_loader=dataloader_mt_t, val_loader=dataloader_mt_v, optimizer=optimizer) # pt.preload_model = True


## META-TEST ##############################################
# create datasets
ds_mt_r_path = os.path.join(opt.data_dir, 'dataloader_mt_r_{}.gz'.format(opt.data_scat.replace('/','_')))
if os.path.exists(ds_mt_r_path):
    dataset_mt_r = compress_pickle.load(ds_mt_r_path, compression="lzma", set_default_extension=False) #gzip
else:
    dataset_mt_r = Data2D(opt, transform=transform_dict['B'], x_files=x_files_mt_r)
    compress_pickle.dump(dataset_mt_r, ds_mt_r_path, compression="lzma", set_default_extension=False) #gzip

# create dataloaders
dataloader_mt_r = DataLoader(dataset=dataset_mt_r,
                                batch_size=len(dataset_mt_r), shuffle=False, drop_last=False,
                                num_workers=1)

for idx, (inp, target, i, xdir, xfn) in enumerate(dataloader_mt_r):
    break

model.eval()

(H, W, C) = (opt.dim, opt.dim, len(opt.x_2D))
img_metas = [{
    'img_shape': (H, W, C),
    'ori_shape': (H, W, C),
    'pad_shape': (H, W, C),
    'filename': xfn__,
    'scale_factor': 1.0,
    'flip': False,
} for xfn__ in xfn]
scores = model.forward(inp, img_metas, gt_semantic_seg=target, return_loss=True)
acc.append([xdir_, float(scores['decode.acc_seg'])])

res = model.inference(inp, img_metas, rescale=False)
val_acc, val_loss, val_losses = validate(val_loader=dataloader_mt_r, model=model, opt=opt)

