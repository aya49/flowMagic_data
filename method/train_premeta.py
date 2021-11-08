import os
import sys
import time
import numpy as np
import torch
import torch.nn as nn
import torch.backends.cudnn as cudnn
import segmentation_models_pytorch as smp

import gc

import tensorboard_logger as tb_logger

from lovasz_losses import lovasz_softmax, iou
from Diceloss import GDiceLossV2 as dice_loss
from dataset import tensor2D3D_

from util import save_checkpoint, load_checkpoint, AverageMeter, adjust_learning_rate
from models import metafreeze_model

def less0_classes(classes, max_class):
    if max_class>1 and classes=='less0':
        return list(range(1,max_class+1))
    elif classes=='present' or classes=='less0':
        return list(range(0,max_class+1))
    return classes

# one epoch validate
def valid_epoch(epoch, val_loader, model, opt, lossfunc, accmetric, classes='present', verbose=True, verboselast=True):
    """One epoch validation"""
    batch_time = AverageMeter()
    losses = AverageMeter()
    top1 = AverageMeter()
    
    ldl = ldl = len(val_loader.dataset.x_files[0])
    with torch.no_grad():
        end = time.time()
        for idx, stuff in enumerate(val_loader):
            # if opt.model == 'setr':
            #     (inp, target, i, xdir, xfn) = stuff
            #     
            #     (H, W, C) = (opt.dim, opt.dim, len(opt.x_2D))
            #     img_metas = [{
            #         'img_shape': (H, W, C),
            #         'ori_shape': (H, W, C),
            #         'pad_shape': (H, W, C),
            #         'filename': xfn_,
            #         'scale_factor': 1.0,
            #         'flip': False,
            #     } for xfn_ in xfn]
            # else:
            (inp, target) = stuff
            
            if torch.cuda.is_available():
                inp = inp.cuda()
                target = target.cuda()
            
            max_class = int(target.max())
            classes = less0_classes(classes, max_class)
            
            # =================== inference =====================
            # if opt.model == 'setr':
            #     scores = model.forward(inp, img_metas, gt_semantic_seg=target, return_loss=True)
            #     # assert isinstance(scores, dict)
            #     acc1 = float(scores['decode.acc_seg'])
            #     loss = float(scores['decode.loss_lovasz'])
            # else:
            
            output = model(inp)
            output = output['out'] if opt.model == 'deeplab3' else output
            # loss = lossfunc(output, target, classes=classes)
            lossfunc = dice_loss()
            # target = tensor2D3D_(target,6).cuda()
            loss = lossfunc.forward(output[:,classes], target)
            acc1 = accmetric(output[:,classes], target)
            
            losses.update(float(loss), inp.size(0))
            top1.update(float(acc1), inp.size(0))
            
            del(inp)
            del(target)
            del(output)
            torch.cuda.empty_cache()
            # gc.collect()
            
            # measure elapsed time
            batch_time.update(time.time() - end)
            end = time.time()
            
            # if idx == ldl-1 and verbose:
            printoutput = verbose and idx==1 if verboselast else verbose
            if printoutput:
                print('Valid [{0}][{1}/{2}]\t'
                      'loss {loss:.3f} ({lossa:.3f})\t'
                      'acc {acc1:.3f} ({acca:.3f})\t'
                      'time {batch_time.val:.2f} ({batch_time.avg:.2f})'.format(
                       epoch, idx, ldl-1, batch_time=batch_time, loss=loss, lossa=losses.avg,
                       acc1=acc1, acca=top1.avg))
    
    return top1.avg, losses.avg


# One epoch training
def train_epoch(epoch, train_loader, model, opt, optimizer, lossfunc, accmetric, classes='present', verbose=True, verboselast=True):
    
    batch_time = AverageMeter()
    data_time = AverageMeter()
    losses = AverageMeter()
    top1 = AverageMeter()
    
    end = time.time()
    ldl = len(train_loader.dataset.x_files[0])
    for idx, stuff in enumerate(train_loader):
        data_time.update(time.time() - end)
        
        # if opt.model == 'setr':
        #     (inp, target, i, xdir, xfn) = stuff
        #     
        #     (H, W, C) = (opt.dim, opt.dim, len(opt.x_2D))
        #     img_metas = [{
        #         'img_shape': (H, W, C),
        #         'ori_shape': (H, W, C),
        #         'pad_shape': (H, W, C),
        #         'filename': xfn_,
        #         'scale_factor': 1.0,
        #         'flip': False,
        #     } for xfn_ in xfn]
        # else:
        (inp, target) = stuff
        
        max_class = int(target.max())
        classes = less0_classes(classes, max_class)
        
        if torch.cuda.is_available():
            inp = inp.cuda()
            target = target.cuda()
        
        # =================== forward =====================
        # if opt.model == 'setr':
        #     ls = model.forward(inp, img_metas, gt_semantic_seg=target, return_loss=True)
        #     loss = ls['decode.loss_lovasz']
        #     acc1 = ls['decode.acc_seg']
        # else:
        output = model(inp) 
        output = output['out'] if opt.model == 'deeplab3' else output
        # loss = lossfunc(output, target, classes=classes)
        lossfunc = dice_loss()
        # target = tensor2D3D_(target,6).cuda()
        loss = lossfunc.forward(output[:,classes], target)
        acc1 = accmetric(output[:,classes], target)
        
        losses.update(float(loss), inp.size(0))
        top1.update(float(acc1), inp.size(0))
        
        # =================== backward =====================
        optimizer.zero_grad()
        loss.backward()
        
        del(inp)
        del(target)
        del(output)
        torch.cuda.empty_cache()
        # gc.collect()
        
        optimizer.step()
        
        # =================== verbose =====================
        batch_time.update(time.time() - end)
        end = time.time()
        
        # print info
        # if idx == ldl-1 and verbose:
        printoutput = verbose and idx==1 if verboselast else verbose
        if printoutput:
            print('Epoch [{0}][{1}/{2}]\t'
                  'loss {loss:.4f} ({lossa:.3f})\t'
                  'acc {acc1:.3f} ({acca:.3f})\t'
                  'time {batch_time.val:.2f} ({batch_time.avg:.2f})'.format(
                  epoch, idx, ldl-1, batch_time=batch_time,
                  data_time=data_time, # 'load {data_time.val:.2f}'
                  acc1=acc1, acca=top1.avg, loss=loss, lossa=losses.avg))
            sys.stdout.flush()
    
    return top1.avg, losses.avg


def train(opt, model, train_loader, val_loader, model_t=None, 
          optimizer=None, lossfunc=None, accmetric=None,
          overwrite=True, classes='less0'):
    
    start = time.time()
    
    optimizer = torch.optim.Adam(model.parameters(), lr=opt.learning_rate, weight_decay=0.0005) if optimizer==None else optimizer
    lossfunc = lovasz_softmax if lossfunc==None else lossfunc
    accmetric = smp.utils.metrics.IoU(threshold=0.5) if accmetric==None else accmetric
    
    # if opt.model == 'setr':
    # else:
    #     optimizer = torch.optim.Adam([dict(params=model.parameters(), lr=opt.learning_rate, weight_decay=0.0005),])
    #     lossfunc = smp.losses.LovaszLoss('multiclass') 
    #     # lossfunc = smp.utils.losses.DiceLoss('multiclass')
    #     # lossfunc = smp.losses.JaccardLoss(mode='multiclass', from_logits=False)
    #     accmetric = [smp.utils.metrics.IoU(threshold=0.5),]
    #     
    #     if opt.mode == 'meta':
    #         model = metafreeze_model(model, opt)
    #         
    #     train_epoch_ = smp.utils.train.TrainEpoch(
    #         model, 
    #         loss=lossfunc, 
    #         metrics=accmetric, 
    #         optimizer=optimizer,
    #         device=torch.device('cuda:0' if opt.gpu else 'cpu'),
    #         verbose=True,
    #     )
    #     valid_epoch_ = smp.utils.train.ValidEpoch(
    #         model, 
    #         loss=lossfunc, 
    #         metrics=accmetric, 
    #         device=torch.device('cuda:0' if opt.gpu else 'cpu'),
    #         verbose=True,
    #     )
    
    if torch.cuda.is_available():
        torch.backends.cudnn.benchmark = True
        if opt.n_gpu > 1:
            model = nn.DataParallel(model)
        # model = model.cuda()
    
    # initialize tensorboard
    logger = tb_logger.Logger(logdir=opt.tb_folder, flush_secs=2)
    
    # save initial checkpoint and load if not overwrite
    epoch_ = 0
    save_file = os.path.join(opt.model_folder, 'ckpt_epoch_{epoch}.pth'.format(epoch=str(epoch_).zfill(3)))
    save_checkpoint(model=model, optimizer=optimizer, save_path=save_file, epoch=epoch_)
    
    if not overwrite:
        ckpts = [x for x in os.listdir(opt.model_folder) if 'ckpt' in x]
        ckpts.sort()
        if '000' not in ckpts[-1]:
            model, _, epoch_ = load_checkpoint(model, os.path.join(opt.model_folder, ckpts[-1]))
        loss_ = np.loadtxt(os.path.join(opt.model_folder, 'loss.csv'), delimiter=',').tolist()
        loss = np.loadtxt(os.path.join(opt.model_folder, 'loss_train.csv'), delimiter=',').tolist()
        acc_ = np.loadtxt(os.path.join(opt.model_folder, 'acc.csv'), delimiter=',').tolist()
        acc = np.loadtxt(os.path.join(opt.model_folder, 'acc_train.csv'), delimiter=',').tolist()
    else:
        acc_, acc, loss_, loss = [], [], [], []
    
    # train: for each epoch
    for epoch in range(epoch_ + 1, opt.epochs + 1):
        
        adjust_learning_rate(epoch, opt, optimizer)
        
        model.train()
        if opt.mode == 'meta':
            model = metafreeze_model(model, opt)
        
        # train!
        time1 = time.time()
        train_acc, train_loss = train_epoch(epoch=epoch, train_loader=train_loader, 
                                            model=model, opt=opt, 
                                            optimizer=optimizer, classes=classes,
                                            lossfunc=lossfunc, accmetric=accmetric,
                                            verbose=epoch%opt.print_freq==0)
        # else:
        #     train_logs  = train_epoch_.run(train_loader)
        #     train_acc = train_logs['dice_loss']
        #     train_loss = train_logs['iou_score']
        time2 = time.time()
        
        logger.log_value('train_acc', train_acc, epoch)
        logger.log_value('train_loss', train_loss, epoch)
        if len(acc_) == 0:
            acc = [train_acc]
            loss = [train_loss]
        else:
            acc.extend([train_acc])
            loss.extend([train_loss])
        
        # validate
        # print('epoch {}, total time {:.2f}'.format(epoch, time2 - time1))
        # regular saving
        if epoch % opt.save_freq == 0:
            # if opt.model == 'setr':
            model.eval()
            val_acc, val_loss = valid_epoch(epoch=epoch, val_loader=val_loader, 
                                            model=model, opt=opt, classes=classes,
                                            lossfunc=lossfunc, accmetric=accmetric,
                                            verbose=epoch%opt.print_freq==0)
            logger.log_value('test_loss', val_loss, epoch)
            logger.log_value('test_acc', val_acc, epoch)
            # else:
            #     valid_logs = valid_epoch_.run(val_loader)
            #     val_acc = valid_logs['dice_loss']
            #     val_loss = valid_logs['iou_score']
            
            if len(acc_) == 0:
                acc_ = [val_acc]
                loss_ = [val_loss]
            else:
                acc_.extend([val_acc])
                loss_.extend([val_loss])
            
            print('==> Saving... {}'.format(opt.model_folder))
            save_file = os.path.join(opt.model_folder, 'ckpt_epoch_{epoch}.pth'.format(epoch=str(epoch).zfill(3)))
            save_checkpoint(model, optimizer, save_file, epoch, opt.n_gpu)
            
            loss_file = os.path.join(opt.model_folder, 'loss.csv')
            np.savetxt(loss_file, loss_, delimiter=', ', fmt="% s")
            loss_file = os.path.join(opt.model_folder, 'loss_train.csv')
            np.savetxt(loss_file, loss, delimiter=', ', fmt="% s")
            
            acc_file = os.path.join(opt.model_folder, 'acc.csv')
            np.savetxt(acc_file, acc_, delimiter=', ', fmt="% s")
            acc_file = os.path.join(opt.model_folder, 'acc_train.csv')
            np.savetxt(acc_file, acc, delimiter=', ', fmt="% s")
        
    
    # save the last model
    save_file = os.path.join(opt.model_folder, '{}_last.pth'.format(opt.model))
    save_checkpoint(model, optimizer, save_file, opt.epochs, opt.n_gpu)
    
    loss_file = os.path.join(opt.model_folder, 'loss.csv')
    np.savetxt(loss_file, loss_, delimiter=', ', fmt="% s")
    loss_file = os.path.join(opt.model_folder, 'loss_train.csv')
    np.savetxt(loss_file, loss, delimiter=', ', fmt="% s")
    
    acc_file = os.path.join(opt.model_folder, 'acc.csv')
    np.savetxt(acc_file, acc_, delimiter=', ', fmt="% s")
    acc_file = os.path.join(opt.model_folder, 'acc_train.csv')
    np.savetxt(acc_file, acc, delimiter=', ', fmt="% s")
    
    dur = time.time() - start
    with open(os.path.join(opt.model_folder, 'time.txt'), mode='a') as file:
        file.write('training/saving time (seconds): %s' % str(dur))
    
    return acc_, loss_, model # i return the model because i like physically seeing it