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

## number of parameters
# sum(p.numel() for p in model.parameters())
## If you want to calculate only the trainable parameters:
# sum(p.numel() for p in model.parameters() if p.requires_grad)

import tensorboard_logger as tb_logger

from GPUtil import showUtilization as gpu_usage # gpu_usage()

from models import create_model

from train import train

print("cuda available")
print(torch.cuda.is_available())

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
x_files_tr = flatx([flatx([[os.path.join(x_den, f) for f in os.listdir(x_den)] for x_den in x_dirs_tr])])
x_files_tr = [x for x in x_files_tr if '__MACOSX' not in x]

x_files_tr_v_ind = random.sample(range(0,len(x_files_tr)), int(len(x_files_tr)/20))
x_files_tr_v = [x_files_tr[x] for x in range(0,len(x_files_tr)) if x in x_files_tr_v_ind]
x_files_tr_t = [x_files_tr[x] for x in range(0,len(x_files_tr)) if x not in x_files_tr_v_ind]

# set some parameters --- if not enough gpu memory, reduce batch_size
opt.num_workers = 32
opt.batch_size = 32
opt.preload_data = True # we pre-load everything so it's faster but takes up more memory
opt.cuda = 'cuda:0'

# create datasets
ds_tr_t_path = os.path.join(opt.data_dir, 'data_tr_t_{}.gz'.format(ds_mt))
if os.path.exists(ds_tr_t_path):
    dataset_tr_t = compress_pickle.load(ds_tr_t_path, compression="lzma", set_default_extension=False) #gzip
else:
    dataset_tr_t = Data2D(opt, transform=transform_dict['A'], x_files=x_files_tr_t)
    compress_pickle.dump(dataset_tr_t, ds_tr_t_path, compression="lzma", set_default_extension=False) #gzip

ds_tr_v_path = os.path.join(opt.data_dir, 'data_tr_v_{}.gz'.format(ds_mt))
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

# train
train(opt, model, dataloader_tr_t, dataloader_tr_v) # opt.preload_model = True


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
## for n_shot in [1, 2, 3, 4, 5, 10, 15, 20]:
n_shot = 10

opt.model_name_meta = '{}_METAdatascat:{}_METAshots:{}'.format(opt.model_name, opt.data_scat, opt.n_shots)

# shot_dir = ./data/x_2Ddensity_euclidean_rankkmed / pregnancy/07_FoxP3CD25_CD4Tcell / 1-5, 10, 15, 20
opt.shot_dir = os.path.join(opt.shot_dir, opt.data_scat, str(opt.n_shots))


                ## META TRAIN ################################################
                # load model
model = create_model(opt).cuda()
ckpt = torch.load(opt.model_path)
model.load_state_dict(torch.load(ckpt)['model'])

# get n-shot samples
mt_dir = os.path.dirname(x_dirs_mts_files[0])
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


