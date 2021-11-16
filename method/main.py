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

# preload data
if opt.preload_data:
    ds_files = []
    for x_dir_mt in x_dirs:#[x for x in x_dirs if 'pregnancy' in x]:
        xdmsplit = x_dir_mt.split('/')
        opt.data_scat = '/'.join(xdmsplit[-2:])
        ds_mt_r_path = os.path.join(opt.data_folder, 'dataloader_mt_r_{}.gz'.format(opt.data_scat.replace('/','_')))
        ds_files.append(ds_mt_r_path)
        if not os.path.exists(ds_mt_r_path):
            print(x_dir_mt)
            x_files_mt = yegz(nomac( [os.path.join(x_dir_mt, f) for f in os.listdir(x_dir_mt)] ))
            dataset_mt_r = Data2D(opt, transform=transform_dict['A'], x_files=x_files_mt)
            compress_pickle.dump(dataset_mt_r, ds_mt_r_path, compression="lzma", set_default_extension=False) #gzip
        # else:
        #     dataset_mt_r = compress_pickle.load(ds_mt_r_path, compression="lzma", set_default_extension=False)
        #     print(x_dir_mt)
        #     print(len(dataset_mt_r))

x_dirs = nomac( flatx([[os.path.join(opt.data_folder, opt.x_2D[0], ds, sc) for 
         sc in os.listdir(os.path.join(opt.data_folder, opt.x_2D[0], ds))] for ds in dss]) )
x_dirs.sort()

## PRE-TRAIN ALL SEQ #############################################
opt.mode = 'pretrain'
baseline = False
basemeta = False
n_shots_baseline = 10

pretrainmode = not baseline
pretrain_all = [[1,2,3],[0,2,3],[0,1,3], [0,1,2]] # if not baseline
meta_all = [[0],[1],[2],[3]] # if not baseline
n_shots = [10,15,5,20,1]

ymask = True

overwrite_pretrain = True
mf = opt.model_folder
ds_files_tr = [x for x in ds_files]
ds_files_tr.sort()

epochs_sample = 500
epochs_pretrain = 100
for ii in range(len(ds_files_tr) if baseline else len(pretrain_all)-1): 
    opt.mode = 'pretrain'
    ds_tr = ''
    if pretrainmode:
        pretrain = pretrain_all[ii]
        meta = meta_all[ii]
        ds_tr = [x for i, x in enumerate(dss) if i in pretrain]
        ds_mt = [x for i, x in enumerate(dss) if i in meta]
    else:
        dscat = ds_files_tr[ii].split('/')[-1].replace('.gz','').replace('dataloader_mt_r_','').replace('_','/',1)
    
    opt.model_folder = '{}:{}'.format(
                        mf.replace(opt.model, '{}{}{}DICE{}{}'.format(
                            opt.model,
                            'BASE' if baseline else 'PRETRAIN',
                            'mask' if ymask else '',
                            '-{}'.format('-'.join(ds_tr) if pretrainmode else ''),
                            '-{}'.format(n_shots_baseline if basemeta else '') if baseline else '')),
                        '{}_{}'.format(str(ii).zfill(2), dscat.replace('/','_')) if baseline else '')
    print('{}: {}'.format(str(ii).zfill(2), opt.model_folder))
    os.makedirs(opt.model_folder, exist_ok=True)
    
    ## initialize model ####
    if 'model' not in locals(): #if i == 0:
        model = create_model(opt).cuda()
        # sum(p.numel() for p in model.parameters())
        model_state = model.state_dict() if opt.n_gpu <= 1 else model.module.state_dict()
    elif baseline:
        model.load_state_dict(model_state)
    
    model_path = os.path.join(opt.model_folder, '{}_last.pth'.format(opt.model))
    if not overwrite_pretrain and os.path.exists(model_path) and pretrainmode:
        model, _, epoch_ = load_checkpoint(model, model_path)
    else:
        # load data
        if baseline:
            if basemeta:
                xdmsplit = x_dirs[ii].split('/')
                opt.data_scat = '/'.join(xdmsplit[-2:])
                x_files_mt = yegz(nomac( [os.path.join(x_dirs[ii], f) for f in os.listdir(x_dirs[ii])] ))
                
                # get n-shot samples
                opt.n_shots = n_shots_baseline
                shot_folder = os.path.join(opt.root_dir, opt.shot_dir, opt.data_scat, str(opt.n_shots))
                x_files_mt_t_ = os.listdir(shot_folder)
                x_files_mt_t = flatx([[x for x in x_files_mt if x_ in x] for x_ in x_files_mt_t_])
                
                dataset_tr_t = Data2D(opt, transform=transform_dict['A'], x_files=x_files_mt_t*(100//len(x_files_mt_t)))
                dataset_tr_v = compress_pickle.load(ds_files_tr[ii], compression="lzma", set_default_extension=False)
            else:
                dataset_tr_t = compress_pickle.load(ds_files_tr[ii], compression="lzma", set_default_extension=False)
                dataset_tr_v = subset_Data2D(dataset_tr_t, len(dataset_tr_t)//10)
        elif pretrainmode:
            ## create datasets ####
            if opt.preload_data:
                ds_files_tr = [x for x in ds_files if any(dsi in x for dsi in [dss[j] for j in pretrain])]
                dataset_tr_t = compress_pickle.load(ds_files_tr[0], compression="lzma", set_default_extension=False) #gzip
                dataset_tr_t.factorize_labels()
                for i in range(1, len(ds_files_tr)):
                    print(ds_files_tr[i])
                    dataset_tr_t_ = compress_pickle.load(ds_files_tr[i], compression="lzma", set_default_extension=False)
                    dataset_tr_t = merge_Data2D( dataset_tr_t, dataset_tr_t_ ) #gzip
                
                dataset_tr_t.factorize_labels()
                dataset_tr_t.transform = transform_dict['A']
            else:
                # train/metatrain data sets denscats folder paths
                x_dirs_tr = nomac( flatx([[os.path.join(opt.data_folder, opt.x_2D[0], ds, sc) for 
                            sc in os.listdir(os.path.join(opt.data_folder, opt.x_2D[0], ds))] for ds in ds_tr]) )
                x_files_tr = yegz(nomac( flatx([flatx([[os.path.join(x_den, f) for f in os.listdir(x_den)] for x_den in x_dirs_tr])]) ))
                
                # split pre-train data set into train (95%) and validation (5%)
                dataset_tr_t = Data2D(opt, transform=transform_dict['A'], x_files=x_files_tr)
            
            dataset_tr_v = subset_Data2D(dataset_tr_t, len(dataset_tr_t)//20)
            dataset_tr_v.transform = transform_dict['B']
        
        # dataset_tr_t.ybig = True
        # dataset_tr_t.ysqueeze = False
        # dataset_tr_v.ysqueeze = False
        tl = len(dataset_tr_t)
        dataset_tr_v.transform = transform_dict['B']
        if hasattr(dataset_tr_t, 'ymask'):
            dataset_tr_t.ymask = ymask
            dataset_tr_v.ymask = ymask
        
        if opt.model == 'deeplab3':
            dataset_tr_t.xcontourdens = True
            dataset_tr_v.xcontourdens = True
        
        dataset_tr_t.loadxy = True
        dataset_tr_v.loadxy = True
        dataset_tr_t.normx = True
        dataset_tr_v.normx = True
        dataloader_tr_t = DataLoader(dataset=dataset_tr_t, sampler=ids(dataset_tr_t), 
                                    batch_size=opt.batch_size, drop_last=True, #shuffle=True, 
                                    num_workers=opt.num_workers)
        dataloader_tr_v = DataLoader(dataset=dataset_tr_v,
                                    batch_size=opt.batch_size, drop_last=False, shuffle=False,
                                    num_workers=opt.num_workers)
        
        # get epochs
        opt.epochs = epochs_sample if baseline else epochs_pretrain
        opt.save_freq = opt.epochs//10
        opt.print_freq = 1
        # opt = update_opt(opt)
        
        # train and validate        
        acc, loss, model = train(opt=opt, model=model, train_loader=dataloader_tr_t, val_loader=dataloader_tr_v, classes='present', overwrite=True) # pt.preload_model = True
        # for par in model.parameters():
        #     print(par)
    
    opt.mode = 'meta'
    mff = opt.model_folder
    if baseline and basemeta:
        n_shots_ = [n_shots_baseline]
    elif baseline:
        n_shots_ = [0]
    elif pretrainmode:
        n_shots_ = n_shots
        x_dirs_mt = nomac( flatx([[os.path.join(opt.data_folder, opt.x_2D[0], ds, sc) for 
                    sc in os.listdir(os.path.join(opt.data_folder, opt.x_2D[0], ds))] for ds in ds_mt]) )
    for n_shot in n_shots_:
        opt.n_shots = n_shot
        for x_dir_mt in [ds_files_tr[ii]] if baseline else x_dirs_mt:
            if pretrainmode:
                xdmsplit = x_dir_mt.split('/')
                opt.data_scat = '/'.join(xdmsplit[-2:])
                x_files_mt = yegz(nomac( [os.path.join(x_dir_mt, f) for f in os.listdir(x_dir_mt)] ))
                
                # get n-shot samples
                opt.model_name_meta = '{}_{}_METAshots:{}'.format(os.path.split(mff)[-1], '_'.join(xdmsplit[-2:]), opt.n_shots) # data_scat e.g. 'pregnancy/07_FoxP3CD25_CD4Tcell'
                shot_folder = os.path.join(opt.root_dir, opt.shot_dir, opt.data_scat, str(opt.n_shots))
                x_files_mt_t_ = os.listdir(shot_folder)
                x_files_mt_t = flatx([[x for x in x_files_mt if x_ in x] for x_ in x_files_mt_t_])
                # x_files_mt_t =  random.sample(x_files_mt, opt.n_shots) ## TEMP!!!!
                # 
                # # get test samples
                # x_files_mt_r = list(set(x_files_mt) - set(x_files_mt_t))
                
                ## META-TRAIN #################################################            
                # create datasets
                dataset_mt_t = Data2D(opt, transform=transform_dict['A'], x_files=x_files_mt_t*(100//len(x_files_mt_t)))
                dataset_mt_v = copy.deepcopy(dataset_mt_t)
                dataset_mt_v.transform = transform_dict['B']
                # if opt.model == 'setr':
                #     dataset_mt_t.loadxy = False
                #     dataset_mt_v.loadxy = False
                # else:
                dataset_mt_t.loadxy = True
                dataset_mt_v.loadxy = True
                if hasattr(dataset_mt_t, 'ymask'):
                    dataset_mt_t.ymask = ymask
                    dataset_mt_v.ymask = ymask
                
                # create dataloaders
                dataloader_mt_t = DataLoader(dataset=dataset_mt_t,# sampler=ids(dataset_mt_t), 
                                            batch_size=min(len(dataset_mt_t.x_files[0]), opt.batch_size), drop_last=True, # shuffle=True, 
                                            num_workers=opt.num_workers)
                dataloader_mt_v = DataLoader(dataset=dataset_mt_v,# sampler=ids(dataset_mt_v), 
                                            batch_size=min(len(dataset_mt_v.x_files[0]), opt.batch_size), drop_last=False, shuffle=False, 
                                            num_workers=opt.num_workers)
                
                # load model
                if 'model' not in locals():
                    model = create_model(opt).cuda()
                ckpt = torch.load(os.path.join(mff, '{}_last.pth'.format(opt.model)))
                # ckpt = torch.load(os.path.join(mff, 'ckpt_epoch_700.pth'))
                model.load_state_dict(ckpt['model'])
                
                # train and validate
                opt.epochs = epochs_sample if baseline else epochs_pretrain
                opt.save_freq = opt.epochs//10
                opt.print_freq = 1
                opt = update_opt(opt)
                opt.model_folder = os.path.join(opt.root_dir, opt.model_dir, opt.model_name_meta)
                os.makedirs(opt.model_folder, exist_ok=True)
                acc, loss, model = train(opt=opt, model=model, train_loader=dataloader_mt_t, val_loader=dataloader_mt_v, classes='present') # pt.preload_model = True
                # acc_path = os.path.join(opt.model_folder, 'acc.csv')
                # loss_path = os.path.join(opt.model_folder, 'loss.csv')
            
            ## META-TEST ##############################################
            # load datasets
            ds_mt_r_path = os.path.join(opt.data_folder, 'dataloader_mt_r_{}.gz'.format(opt.data_scat.replace('/','_')))
            
            if baseline and basemeta:
                dataset_mt_r = dataset_tr_v
            elif baseline and not basemeta:
                dataset_mt_r = dataset_tr_t
            elif pretrainmode:
                dataset_mt_r = compress_pickle.load(ds_mt_r_path, compression="lzma", set_default_extension=False) #gzip
            
            dataset_mt_r.transform = transform_dict['B']
            dataset_mt_r.loadxy = False
            if hasattr(dataset_mt_r, "ymask"):
                dataset_mt_r.ymask = ymask
            
            # create dataloaders
            dataloader_mt_r = DataLoader(dataset=dataset_mt_r,
                                batch_size=10, shuffle=False, drop_last=False,
                                num_workers=opt.num_workers)
            
            model.eval()
            total_r = len(dataset_mt_r)
            res_dir = os.path.join(opt.data_folder.replace('/data/','/results/'), 'method/{}/{}/{}'.format( os.path.split(mff)[-1].split('_')[0], opt.n_shots, opt.data_scat))
            os.makedirs(res_dir, exist_ok=True)
            
            for idx, stuff in enumerate(dataloader_mt_r):
                (inp, target, _, xdir, xfn) = stuff
                
                # if opt.model == 'setr':
                #     inp, target, img_metas = prep_input(inp, target, xfn)
                #     res = model.inference(inp, img_metas, rescale=False)
                # else:
                if torch.cuda.is_available():
                    inp = inp.cuda()
                    # target = target.cuda()
                if opt.model == 'setr':
                    res = model(inp)
                elif opt.model == 'deeplab3':
                    res = model(inp)['out']
                else:
                    res = model.predict(inp)
                
                for xfi in range(len(xfn)):
                    res_ = res[xfi].squeeze()
                    res_vals, res_ind = torch.max(res_, 0) # 3D to 2D
                    res_ind[inp[xfi][0].squeeze()==0] = 0
                    res_ind = pd.DataFrame(res_ind.cpu().detach().numpy())
                    
                    res_file = os.path.join(res_dir, xfn[xfi]) # ends with gz so auto compress
                    res_ind.to_csv(res_file, index=False, header=False, compression='gzip')
                
                del(inp)
                del(res)
                torch.cuda.empty_cache()

