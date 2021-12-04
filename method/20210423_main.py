## REFERENCES #############################################
# for image classification: https://github.com/WangYueFt/rfs
# transunet: https://github.com/Beckschen/TransUNet/tree/main/networks


## MODULES ################################################
# this training was done on 1 gpu and 32 workers
'''
conda activate pytorch_seg

module load LANG/PYTHON/3.7.6
module load LIB/OPENCV/3.4.9-PY376-CUDA
module load TOOLS/PYTORCH/1.7.1-CUDA101-PY376 # automatically loads CUDA 10.1
nvcc --version # default CUDA 9.1 # CUDA version

cd flowMagic_data/src/method
'''

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

import cv2

# torch
import torch
from torch.utils.data import DataLoader
import torchvision.transforms as tr
# from torchviz import make_dot # creates an image of the model
from torchsampler import ImbalancedDatasetSampler as ids # pip install https://github.com/ufoym/imbalanced-dataset-sampler/archive/master.zip

from opt import parse_options, update_opt
from util import prep_input, visualize, load_checkpoint, nomac, yegz, flatx
from transform import transform_dict
from dataset import Data2D, merge_Data2D, subset_Data2D
from models import create_model
from train_premeta import train

from Diceloss import GDiceLossV2 as dice_loss
from Diceloss import BinaryDiceLoss as dice_loss_binary

## number of model parameters
# sum(p.numel() for p in model.parameters())
## If you want to calculate only the trainable parameters:
# sum(p.numel() for p in model.parameters() if p.requires_grad)

print("cuda available")
print(torch.cuda.is_available())

# options
opt = parse_options()
mf = opt.model_folder

## DATA: datasets x 4 ########################################

# dataset names ['HIPCmyeloid', 'sangerP2', 'HIPCbcell', 'pregnancy']
dss = nomac( os.listdir(os.path.join(opt.data_folder, opt.x_2D[0])) )

# x dataset directories ['~/data/2D/x_2Dscatter/HIPCmyeloid/CD3SSCA_HLADR-CD14-_', ...]
x_dirs = nomac( flatx([[os.path.join(opt.data_folder, opt.x_2D[0], ds, sc) for 
                sc in os.listdir(os.path.join(opt.data_folder, opt.x_2D[0], ds))] for ds in dss]) )
x_dirs.sort()

# preload data; alternatively run file load_data.py to load everything
if opt.preload_data:
    ds_files = []
    for x_dir in x_dirs: #[x for x in x_dirs if 'pregnancy' in x]:
        data_scat = '/'.join(x_dir.split('/')[-2:]) # 'HIPCmyeloid/CD3SSCA_HLADR-CD14-_'
        ds_data_path = os.path.join(opt.data_folder, 
                                    'dataloader_mt_r_{}.gz'.format(data_scat.replace('/','_'))) # '~/data/2D/dataloader_mt_r_HIPCmyeloid_CD3SSCA_HLADR-CD14-_.gz'
        ds_files.append(ds_data_path)
        if not os.path.exists(ds_data_path):
            print(x_dir)
            x_files = yegz(nomac( [os.path.join(x_dir, f) for f in os.listdir(x_dir)] ))
            dataset = Data2D(opt, transform=transform_dict['A'], x_files=x_files)
            compress_pickle.dump(dataset, ds_data_path, compression="lzma", 
                                 set_default_extension=False) #gzip
        # # test loading of data sets
        # else:
        #     dataset = compress_pickle.load(ds_data_path, compression="lzma", set_default_extension=False)
        #     print(x_dir)
        #     print(len(dataset))

## PARAMETERS #############################################
baseline = True # no pre-training
basemeta = True # if baseline, train with k samples, else train with all samples
n_shots_baseline = [10]
epochs_sample = 500

pretrainmode = not baseline
n_shots = [10,15,5,20,1]
epochs_pretrain = 500
epochs_metatrain = 200
pretrain_all = [[0,1,2], [1,2,3],[0,2,3],[0,1,3]] # if not baseline
meta_all = [[3],[0],[1],[2]] # if not baseline

ymask = True
singlecpop = True
weightbg0 = False # weight background as 0 for calculating loss
v_ratio = 20 # validation set ratio = len(training)//v_ratio
overwrite_pretrain = True

ds_tr = ''
opt.preload_data = True # we pre-load everything so it's faster but takes up more memory
opt.num_workers = 32
opt.batch_size = 32 # if not enough gpu memory, reduce batch_size

for ii in range(len(ds_files) if baseline else len(pretrain_all)-1): #[x for x in range(len(ds_files)) if 'pregnancy' in ds_files[x]]:
    opt.mode = 'pretrain'
    if pretrainmode:
        ds_tr = [x for i, x in enumerate(dss) if i in pretrain_all[ii]]
        ds_mt = [x for i, x in enumerate(dss) if i in meta_all[ii]]
    else:
        dscat = ds_files[ii].split('/')[-1].replace('.gz','').replace('dataloader_mt_r_','').replace('_','/',1)
    
    opt.model_folder = '{}:{}'.format(
                        mf.replace(opt.model, '{}{}{}DICE{}{}{}'.format(
                            opt.model,
                            'BASE' if baseline else 'PRETRAIN',
                            'mask' if ymask else 'raw',
                            'wbg0' if weightbg0 else 'wbg1',
                            'Singlecpop' if singlecpop else 'Multicpop',
                            '-{}'.format('-'.join(ds_tr) if pretrainmode else ''))),
                        '{}_{}'.format(str(ii).zfill(2), dscat.replace('/','_')) if baseline else '')
    print('{}: {}'.format(str(ii).zfill(2), opt.model_folder))
    os.makedirs(opt.model_folder, exist_ok=True)
    
    ## initialize model ####
    if 'model' not in locals():
        model = create_model(opt, singlecpop).cuda() # sum(p.numel() for p in model.parameters())
        model_state = model.state_dict() if opt.n_gpu <= 1 else model.module.state_dict()
    elif baseline:
        model.load_state_dict(model_state)
    
    model_path = os.path.join(opt.model_folder, '{}_last.pth'.format(opt.model))
    if not overwrite_pretrain and os.path.exists(model_path) and pretrainmode:
        model, _, epoch_ = load_checkpoint(model, model_path)
    elif pretrainmode:
        ## create datasets ####
        if opt.preload_data:
            ds_files_tr_ = [x for x in ds_files if any(dsi in x for dsi in [dss[j] for j in pretrain_all[ii]])]
            dataset_tr_t = compress_pickle.load(ds_files_tr_[0], compression="lzma", 
                                                set_default_extension=False) #gzip
            for i in range(1, len(ds_files_tr_)):
                print(ds_files_tr_[i])
                dataset_tr_t_ = compress_pickle.load(ds_files_tr_[i], compression="lzma", 
                                                     set_default_extension=False)
                dataset_tr_t = merge_Data2D( dataset_tr_t, dataset_tr_t_ ) #gzip
            
            dataset_tr_t.factorize_labels()
        else:
            # train/metatrain data sets denscats folder paths
            x_dirs_tr = nomac( flatx([[os.path.join(opt.data_folder, opt.x_2D[0], ds, sc) for 
                               sc in os.listdir(os.path.join(opt.data_folder, opt.x_2D[0], ds))] for 
                               ds in ds_tr]) )
            x_files_tr = yegz(nomac( flatx([flatx([[os.path.join(x_den, f) for 
                                     f in os.listdir(x_den)] for 
                                     x_den in x_dirs_tr])]) ))
            
            dataset_tr_t = Data2D(opt, transform=transform_dict['A'], x_files=x_files_tr)
        
        # dataset_tr_t.ybig = True
        # dataset_tr_t.ysqueeze = False
        # dataset_tr_t.ymask = ymask
        # dataset_tr_t.loadxy = True
        dataset_tr_t.normx = True
        dataset_tr_t.cpop = -1 if singlecpop else None
        dataset_tr_t.transform = transform_dict['A']
        
        # split pre-train data set into train (95%) and validation (5%)
        tl = len(dataset_tr_t)
        if opt.preload_data:
            dataset_tr_v = compress_pickle.load(ds_files_tr_[0], compression="lzma", 
                                                set_default_extension=False) #gzip
        
        dataset_tr_v = subset_Data2D(dataset_tr_t, tl//v_ratio, 
                                     dataset_tr_v if opt.preload_data else None)
        dataset_tr_v.transform = transform_dict['B']
        
        dataloader_tr_t = DataLoader(dataset=dataset_tr_t, sampler=ids(dataset_tr_t), 
                                    batch_size=opt.batch_size, drop_last=True, #shuffle=True, 
                                    num_workers=opt.num_workers)
        dataloader_tr_v = DataLoader(dataset=dataset_tr_v,
                                    batch_size=opt.batch_size, drop_last=False, shuffle=False,
                                    num_workers=opt.num_workers)
        
        ## pre-train model ####
        opt.epochs = epochs_pretrain
        opt.save_freq = max(1, opt.epochs//10)
        opt.print_freq = 1
        acc, loss, model = train(opt=opt, model=model, classes='present', overwrite=True, 
                                 train_loader=dataloader_tr_t, val_loader=dataloader_tr_v,
                                 lossfunc=dice_loss_binary() if singlecpop else dice_loss(),
                                 weightbg0=weightbg0) # pt.preload_model = True
        # for par in model.parameters():
        #     print(par)
    
    opt.mode = 'meta'
    mff = opt.model_folder
    if baseline:
        n_shots_ = n_shots_baseline if basemeta else [0] 
    elif pretrainmode:
        n_shots_ = n_shots
        x_dirs_mt = nomac( flatx([[os.path.join(opt.data_folder, opt.x_2D[0], ds, sc) for 
                    sc in os.listdir(os.path.join(opt.data_folder, opt.x_2D[0], ds))] for 
                    ds in ds_mt]) )
    for n_shot in n_shots_:
        opt.n_shots = n_shot
        for x_dir_mt in [x_dirs[ii]] if baseline else x_dirs_mt:
            xdmsplit = x_dir_mt.split('/')
            opt.data_scat = '/'.join(xdmsplit[-2:])
            
            opt.model_folder = '{}_{}_METAshots:{}'.format(mff, '_'.join(xdmsplit[-2:]) if pretrainmode else '', opt.n_shots)
            os.makedirs(opt.model_folder, exist_ok=True)
            
            ## create training dataset ####
            if baseline and not basemeta:
                dataset_mt_t = compress_pickle.load(ds_files[ii], compression="lzma", 
                                                    set_default_extension=False)
                dataset_mt_v = subset_Data2D(dataset_mt_t, len(dataset_mt_t)//v_ratio)
            else:
                x_files_mt = yegz(nomac( [os.path.join(x_dir_mt, f) for 
                                         f in os.listdir(x_dir_mt)] ))
                
                # get n-shot samples
                shot_folder = os.path.join(opt.root_dir, opt.shot_dir, opt.data_scat, str(opt.n_shots))
                x_files_mt_t_ = os.listdir(shot_folder)
                x_files_mt_t = flatx([[x for x in x_files_mt if x_ in x] for x_ in x_files_mt_t_])
                # x_files_mt_t =  random.sample(x_files_mt, opt.n_shots) ## TEMP!!!!
                # 
                # # get test samples
                # x_files_mt_r = list(set(x_files_mt) - set(x_files_mt_t))
                
                ## create datasets ####
                dataset_mt_t = Data2D(opt, transform=transform_dict['A'], x_files=x_files_mt_t*(100//len(x_files_mt_t)))
                
                dataset_mt_v = copy.deepcopy(dataset_mt_t)
            
            dataset_mt_t.ymask = ymask
            dataset_mt_t.loadxy = True
            dataset_mt_t.normx = True
            dataset_mt_v.transform = transform_dict['B']
            dataset_mt_v.ymask = ymask
            dataset_mt_v.loadxy = True
            dataset_mt_v.normx = True
            
            # train and validate
            opt.epochs = epochs_sample if baseline else epochs_metatrain
            opt.save_freq = max(1, opt.epochs//10)
            opt.print_freq = 1
            
            num_class = int(dataset_mt_t.y[0].max())
            for cpop in range(1,num_class) if singlecpop else [0]:
                if cpop>0:
                    dataset_mt_t.cpop = 0
                    y = dataset_mt_t.__getitem__(0)[1]
                    xind, yind, w, h = cv2.boundingRect(np.uint8(y[0] == cpop))
                    if len(dataset_mt_t)>1:
                        xind2 = xind+w
                        yind2 = yind+h
                        for di in range(len(dataset_mt_t)):
                            x, y = dataset_mt_t.__getitem__(di)
                            xind_, yind_, w_, h_ = cv2.boundingRect(np.uint8(y[0] == cpop))
                            xind = min(xind, xind_)
                            yind = min(yind, yind_)
                            xind2 = max(xind_+w_, xind2)
                            yind2 = max(yind_+h_, yind2)
                        w = xind2-xind
                        h = yind2-yind
                    cpopdim = [xind, yind, w, h]
                    tr_resize = tr.Resize((w, h))
                    
                    dataset_mt_t.dim = cpopdim
                    dataset_mt_v.dim = cpopdim
                
                dataset_mt_t.transform = transform_dict['A']
                
                # create dataloaders
                dataset_mt_t.cpop = cpop
                dataset_mt_v.cpop = cpop
                dataloader_mt_t = DataLoader(dataset=dataset_mt_t,# sampler=ids(dataset_mt_t), 
                                            batch_size=min(len(dataset_mt_t.x_files[0]), opt.batch_size), drop_last=True, # shuffle=True, 
                                            num_workers=opt.num_workers)
                dataloader_mt_v = DataLoader(dataset=dataset_mt_v,# sampler=ids(dataset_mt_v), 
                                            batch_size=min(len(dataset_mt_v.x_files[0]), opt.batch_size), drop_last=False, shuffle=False, 
                                            num_workers=opt.num_workers)
                
                # load model
                if pretrainmode:
                    ckpt = torch.load(os.path.join(mff, '{}_last.pth'.format(opt.model)))
                    # ckpt = torch.load(os.path.join(mff, 'ckpt_epoch_700.pth'))
                    model.load_state_dict(ckpt['model'])
                else:
                    model.load_state_dict(model_state)
                
                acc, loss, model = train(opt=opt, model=model, classes='present', overwrite=True, 
                                        train_loader=dataloader_mt_t, val_loader=dataloader_mt_v,
                                        lossfunc=dice_loss_binary() if singlecpop else dice_loss(),
                                 weightbg0=weightbg0) # pt.preload_model = True
                # for par in model.parameters():
                #     print(par)
                # acc_path = os.path.join(opt.model_folder, 'acc.csv')
                # loss_path = os.path.join(opt.model_folder, 'loss.csv')
                
                ## META-TEST ##############################################
                # load datasets
                ds_mt_r_path = os.path.join(opt.data_folder, 'dataloader_mt_r_{}.gz'.format(opt.data_scat.replace('/','_')))
                
                dataset_mt_r = compress_pickle.load(ds_mt_r_path, compression="lzma", set_default_extension=False) #gzip
                
                dataset_mt_r.transform = transform_dict['B']
                dataset_mt_r.loadxy = False
                dataset_mt_r.ymask = ymask
                
                # create dataloaders
                if cpop>0:
                    dataset_mt_r.dim = cpopdim
                dataset_mt_r.cpop = cpop
                dataset_mt_r.dim = None if cpop==0 else cpopdim
                # dataloader_mt_r = DataLoader(dataset=dataset_mt_r,
                #                           batch_size=10, shuffle=False, drop_last=False,
                #                           num_workers=0)
                
                model.eval()
                total_r = len(dataset_mt_r)
                res_dir = os.path.join(opt.data_folder.replace('/data/','/results/'), 'method/{}/{}/{}'.format( os.path.split(mff)[-1].split('_')[0], opt.n_shots, opt.data_scat))
                os.makedirs(res_dir, exist_ok=True)
                
                endclass = cpop>0 and cpop==num_class
                # for idx, stuff in enumerate(dataloader_mt_r):
                for ids in range(len(dataset_mt_r)):
                    inp_, target, _, xdir, xfn = dataset_mt_r.__getitem__(ids)
                    inp = inp_.unsqueeze(0)
                    
                    # if opt.model == 'setr':
                    #     inp, target, img_metas = prep_input(inp, target, xfn)
                    #     res = model.inference(inp, img_metas, rescale=False)
                    # else:
                    inp = inp.cuda() if torch.cuda.is_available() else inp
                    
                    if opt.model == 'setr':
                        res = model(inp)
                    elif opt.model == 'deeplab3':
                        res = model(inp)['out']
                    else:
                        res = model.predict(inp)
                    
                    res_file = os.path.join(res_dir, xfn) # ends with gz so auto compress
                    if cpop==0:
                        res_ = res.squeeze()
                    else:
                        res_t_ = tr_resize(res).squeeze()
                        res_t = torch.zeros(opt.dim, opt.dim)
                        res_t[xind:(xind+w),yind:(yind+h)] = res_t_
                        if cpop==1:
                            if endclass:
                                res_ind = res_t.round().int() # i just like seeing assignments
                            else:
                                compress_pickle.dump([res_t], '{}_temp.gz'.format(res_file), compression="lzma", set_default_extension=False) #gzip
                        else:
                            res_temp_ = res_t
                            res_temp = compress_pickle.load('{}_temp.gz'.format(res_file), compression="lzma", set_default_extension=False)
                            res_temp.append(res_temp_)
                            compress_pickle.dump(res_temp, '{}_temp.gz'.format(res_file), compression="lzma", set_default_extension=False) #gzip
                            if endclass:
                                res_ = torch.stack(res_temp)
                                max_class, mcind = torch.max(res_, 0)
                                max_class = max_class<.5
                                max_class = max_class.int()
                                res_ = torch.stack([max_class]+res_temp)
                    
                    if cpop==0 or (cpop>1 and endclass):
                        res_vals, res_ind = torch.max(res_, 0) # 3D to 2D
                        res_ind[inp[0][0].squeeze()==0] = 0
                        
                        res_ind = pd.DataFrame(res_ind.cpu().detach().numpy())
                        res_ind.to_csv(res_file, index=False, header=False, compression='gzip')
                    
                    del(inp)
                    del(res)
                    torch.cuda.empty_cache()

