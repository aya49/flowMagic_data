import os
import sys
import time
import torch
import torch.nn as nn
import torch.backends.cudnn as cudnn
import segmentation_models_pytorch as smp

import gc

import tensorboard_logger as tb_logger

from mmseg.models.losses import lovasz_loss as ll
from util import save_checkpoint, load_checkpoint, AverageMeter, adjust_learning_rate
from models import metafreeze_model

def valid_epoch(val_loader, model, opt):
    """One epoch validation"""
    batch_time = AverageMeter()
    losses = AverageMeter()
    top1 = AverageMeter()
    
    with torch.no_grad():
        end = time.time()
        for idx, (inp, target, i, xdir, xfn) in enumerate(val_loader):
            
            inp = inp.float()
            target = target.float()
            if torch.cuda.is_available():
                inp = inp.cuda()
                target = target.cuda()
            
            (H, W, C) = (opt.dim, opt.dim, len(opt.x_2D))
            img_metas = [{
                'img_shape': (H, W, C),
                'ori_shape': (H, W, C),
                'pad_shape': (H, W, C),
                'filename': xfn_,
                'scale_factor': 1.0,
                'flip': False,
            } for xfn_ in xfn]
            
            # compute output
            scores = model.forward(inp, img_metas, gt_semantic_seg=target, return_loss=True)
            # assert isinstance(scores, dict)
            # loss = criterion(output, target)
            
            # measure accuracy and record loss
            # acc1, acc5 = accuracy(output, target, topk=(1, 5))
            acc1 = float(scores['decode.acc_seg'])
            loss = float(scores['decode.loss_lovasz'])
            losses.update(loss, inp.size(0))
            top1.update(acc1, inp.size(0))
            
            # measure elapsed time
            batch_time.update(time.time() - end)
            end = time.time()
            
            if idx % opt.print_freq == 0:
                print('Test: [{0}/{1}]\t'
                      'Time {batch_time.val:.3f} ({batch_time.avg:.3f})\t'
                      'Loss ({loss:.4f})\t'
                      'Acc@1 ({acc1:.3f})'.format(
                       idx, len(val_loader), batch_time=batch_time, loss=loss,
                       acc1=acc1))

        print(' * Acc@1 {acc1:.3f}'.format(acc1=top1.avg))
    
    return top1.avg, losses.avg
    


def train_epoch(epoch, train_loader, model, optimizer, opt):
    # One epoch training
    
    set_cuda = torch.cuda.is_available()
    
    # # set modules as train()
    # if opt.mode == 'distill':
    #     # set teacher as eval()
    #     model_t = model[-1]
    #     model_t.train()
    #     model_t.eval()
    
    #     model = model[0]
    
    #     criterion_cls = criterion[0]
    #     criterion_div = criterion[1]
    #     criterion_kd = criterion[2]
    
    batch_time = AverageMeter()
    data_time = AverageMeter()
    losses = AverageMeter()
    top1 = AverageMeter()
    
    end = time.time()
    for idx, (inp, target, i, xdir, xfn) in enumerate(train_loader):
        data_time.update(time.time() - end)
        
        # if opt.mode == 'distill':
        #     inp, target, idx, _ = enum
        #     if opt.distill in ['contrast']:
        #         inp, target, index, contrast_idx = enum
        #         contrast_idx = contrast_idx.cuda() if set_cuda else contrast_idx
        #     else:
        #         inp, target, idx, _ = data
        #     index = index.cuda() if set_cuda else index
        # else:
        #     inp, target, idx, _ = enum

        inp = inp.float()
        inp = inp.cuda() if set_cuda else inp
        target = target.cuda() if set_cuda else target
        
        (H, W, C) = (opt.dim, opt.dim, len(opt.x_2D))
        img_metas = [{
            'img_shape': (H, W, C),
            'ori_shape': (H, W, C),
            'pad_shape': (H, W, C),
            'filename': xfn_,
            'scale_factor': 1.0,
            'flip': False,
        } for xfn_ in xfn]
        
        # ===================forward=====================
        ls = model.forward(inp, img_metas, gt_semantic_seg=target, return_loss=True)
        loss = float(ls['decode.loss_lovasz']) # ll.lovasz_softmax(output, target, classes='present', per_image=True)
        acc1 = float(ls['decode.acc_seg']) # ll.iou(output, target, opt.n_class, EMPTY=1., ignore=None, per_image=True)
        
        # outputs = model.train_step(dict(img=inp, img_metas=img_metas, gt_semantic_seg=target), None)
        
        losses.update(loss, inp.size(0))
        top1.update(acc1, inp.size(0))
        
        # ===================backward=====================
        optimizer.zero_grad()
        ls['decode.loss_lovasz'].backward()
        
        del(inp)
        del(target)
        torch.cuda.empty_cache()
        # gc.collect()
        
        optimizer.step()
        
        # ===================meters=====================
        batch_time.update(time.time() - end)
        end = time.time()
        
        # tensorboard logger
        pass
        
        # print info
        if idx % opt.print_freq == 0:
            print('Epoch: [{0}][{1}/{2}]\t'
                    'Time {batch_time.val:.3f} ({batch_time.avg:.3f})\t'
                    'Data {data_time.val:.3f} ({data_time.avg:.3f})\t'
                    'Loss ({loss:.4f})\t'
                    'Acc@1 ({acc1:.3f})\t'.format(
                    epoch, idx, len(train_loader), batch_time=batch_time,
                    data_time=data_time, loss=loss, acc1=acc1))
            sys.stdout.flush()
        
    print(' * Acc@1 {acc1:.3f}'.format(acc1=top1.avg))
    
    return top1.avg, losses.avg


def train(opt, model, train_loader, val_loader, optimizer, model_t=None):
    
    if opt.model == 'setr':
        optimizer = torch.optim.Adam(model.parameters(), lr=opt.learning_rate, weight_decay=0.0005)
    else:
        optimizer = torch.optim.Adam([dict(params=model.parameters(), lr=opt.learning_rate, weight_decay=0.0005),])
        loss = smp.losses.LovaszLoss(mode='multiclass') # loss = smp.utils.losses.DiceLoss()
        metrics = [smp.utils.metrics.IoU(threshold=0.5),]
        
        train_epoch_ = smp.utils.train.TrainEpoch(
            model, 
            loss=loss, 
            metrics=metrics, 
            optimizer=optimizer,
            device=torch.device('cuda:0' if opt.gpu else 'cpu'),
            verbose=True,
        )
        valid_epoch_ = smp.utils.train.ValidEpoch(
            model, 
            loss=loss, 
            metrics=metrics, 
            device=torch.device('cuda:0' if opt.gpu else 'cpu'),
            verbose=True,
        )
    
    if torch.cuda.is_available():
        torch.backends.cudnn.benchmark = True
        if opt.n_gpu > 1:
            model = nn.DataParallel(model)
        # model = model.cuda()
        # ll.cuda()

        
    
    # tensorboard
    logger = tb_logger.Logger(logdir=opt.tb_folder, flush_secs=2)
    
    # set cosine annealing scheduler
    if opt.cosine:
        eta_min = opt.learning_rate * (opt.lr_decay_rate ** 3)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, opt.epochs, eta_min, -1)
    
    # train
    epoch_ = 0
    save_file = os.path.join(opt.model_folder, 'ckpt_epoch_{epoch}.pth'.format(epoch=str(epoch_).zfill(3)))
    save_checkpoint(model=model, optimizer=optimizer, save_path=save_file, epoch=epoch_)
    
    ckpts = os.listdir(opt.model_folder)
    ckpts.sort()
    print(ckpts)
    if 'ckpt' in ckpts[-1] and '000' not in ckpts[-1]:
        model, optimizer, epoch_ = load_checkpoint(model, os.path.join(opt.model_folder, ckpts[-1]))
    
    print(epoch_)
    
    acc = []
    loss = []
    for epoch in range(epoch_ + 1, opt.epochs + 1):
        
        if opt.cosine:
            scheduler.step()
        else:
            adjust_learning_rate(epoch, opt, optimizer)
        
        print("==> training")
        time1 = time.time()
        if opt.mode == 'meta':
            model = metafreeze_model(model, opt)
        if opt.model == 'setr':
            model.train()
            train_acc, train_loss = train_epoch(epoch=epoch, train_loader=train_loader, model=model, optimizer=optimizer, opt=opt)
        else:
            train_logs  = train_epoch_.run(train_loader)
        time2 = time.time()
        print('epoch {}, total time {:.2f}'.format(epoch, time2 - time1))
        
        logger.log_value('train_acc', train_acc, epoch)
        logger.log_value('train_loss', train_loss, epoch)
        

        if opt.model == 'setr':
            model.eval()
            val_acc, val_loss = valid_epoch(val_loader=val_loader, model=model, opt=opt)
            acc = acc.append(val_loss)
            loss = loss.append(val_acc)
        else:
            valid_logs = valid_epoch_.run(val_loader)
            acc = acc.append(valid_logs['dice_loss'])
            loss = loss.append(valid_logs['iou_score'])
        logger.log_value('test_acc', val_acc, epoch)
        logger.log_value('test_loss', val_loss, epoch)
        
        # regular saving
        if epoch % opt.save_freq == 0:
            print('==> Saving...')
            save_file = os.path.join(opt.model_folder, 'ckpt_epoch_{epoch}.pth'.format(epoch=str(epoch).zfill(3)))
            save_checkpoint(model, optimizer, save_file, epoch, opt.n_gpu)
    
    # save the last model
    save_file = os.path.join(opt.model_folder, '{}_last.pth'.format(opt.model))
    save_checkpoint(model, optimizer, save_file, opt.epochs, opt.n_gpu)
    
    return acc, loss, model # yes, i return the model because i like seeing it there