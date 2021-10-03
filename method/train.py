import os
import sys
import time
import torch
import torch.nn as nn
import torch.backends.cudnn as cudnn

import tensorboard_logger as tb_logger

from loss import lovasz_softmax as ll

from util import save_checkpoint, load_checkpoint, AverageMeter, adjust_learning_rate


def validate(val_loader, model, criterion, opt, accuracy):
    """One epoch validation"""
    batch_time = AverageMeter()
    losses = AverageMeter()
    top1 = AverageMeter()

    # switch to evaluate mode
    model.eval()

    with torch.no_grad():
        end = time.time()
        for idx, (inp, target, i, xdir, xfn) in enumerate(val_loader):

            inp = inp.float()
            if torch.cuda.is_available():
                inp = inp.cuda()
                target = target.cuda()

            img_metas = [{
                'filename': xfn_
            } for xfn_ in xfn]

            # compute output
            output = model(inp, img_metas)
            loss = criterion(output, target)

            # measure accuracy and record loss
            acc1, acc5 = accuracy(output, target, topk=(1, 5))
            losses.update(loss.item(), inp.size(0))
            top1.update(acc1[0], inp.size(0))

            # measure elapsed time
            batch_time.update(time.time() - end)
            end = time.time()

            if idx % opt.print_freq == 0:
                print('Test: [{0}/{1}]\t'
                      'Time {batch_time.val:.3f} ({batch_time.avg:.3f})\t'
                      'Loss {loss.val:.4f} ({loss.avg:.4f})\t'
                      'Acc@1 {top1.val:.3f} ({top1.avg:.3f})\t'
                      'Acc@5 {top5.val:.3f} ({top5.avg:.3f})'.format(
                       idx, len(val_loader), batch_time=batch_time, loss=losses,
                       top1=top1))

        print(' * Acc@1 {top1.avg:.3f}'.format(top1=top1))

    return top1.avg, losses.avg


def train_epoch(epoch, train_loader, model, criterion, optimizer, opt):
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
    
    model.train()

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
        
        # inp = inp.float()
        inp = inp.cuda() if set_cuda else inp
        target = target.cuda() if set_cuda else target

        img_metas = [{
            'filename': xfn_
        } for xfn_ in xfn]

        # ===================forward=====================
        if opt.mode == 'distill':
            preact = False
            # if opt.distill in ['abound', 'overhaul']:
            #     preact = True
            # feat, logit = model(inp, is_feat=True)
            # with torch.no_grad():
            #     feat_t, logit_t = model_t(inp, is_feat=True)
            #     feat_t = [f.detach() for f in feat_t]

            # # cls + kl div
            # loss_cls = criterion_cls(logit, target)
            # loss_div = criterion_div(logit, logit_t)

            # # other kd beyond KL divergence
            # if opt.distill == 'kd':
            #     loss_kd = 0
            # elif opt.distill == 'contrast':
            #     f = module_list[1](feat[-1])
            #     f_t = module_list[2](feat_t[-1])
            #     loss_kd = criterion_kd(f, f_t, index, contrast_idx)
            # elif opt.distill == 'hint':
            #     f = feat[-1]
            #     f_t = feat_t[-1]
            #     loss_kd = criterion_kd(f, f_t)
            # elif opt.distill == 'attention':
            #     g = feat[1:-1]
            #     g_t = feat_t[1:-1]
            #     loss_group = criterion_kd(g, g_t)
            #     loss_kd = sum(loss_group)
            # else:
            #     raise NotImplementedError(opt.distill)

            # loss = opt.gamma * loss_cls + opt.alpha * loss_div + opt.beta * loss_kd 
            # acc1 = ll.iou(output, target, opt.n_class, EMPTY=1., ignore=None, per_image=True)
        else:
            output = model.forward(inp, img_metas)
            loss = ll.lovasz_softmax(output, target, classes='present', per_image=True, ignore=None)
            acc1 = ll.iou(output, target, opt.n_class, EMPTY=1., ignore=None, per_image=True)
        
        losses.update(loss.item(), inp.size(0))
        top1.update(acc1[0], inp.size(0))

        # ===================backward=====================
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        # ===================meters=====================
        batch_time.update(time.time() - end)
        end = time.time()

        # print info
        if idx % opt.print_freq == 0:
            print('Epoch: [{0}][{1}/{2}]\t'
                    'Time {batch_time.val:.3f} ({batch_time.avg:.3f})\t'
                    'Data {data_time.val:.3f} ({data_time.avg:.3f})\t'
                    'Loss {loss.val:.4f} ({loss.avg:.4f})\t'
                    'Acc@1 {top1.val:.3f} ({top1.avg:.3f})\t'
                    'Acc@5 {top5.val:.3f} ({top5.avg:.3f})'.format(
                    epoch, idx, len(train_loader), batch_time=batch_time,
                    data_time=data_time, loss=losses, top1=top1))
            sys.stdout.flush()

    print(' * Acc@1 {top1.avg:.3f} Acc@5 {top5.avg:.3f}'.format(top1=top1))

    return top1.avg, losses.avg


def train(opt, model, train_loader, val_loader, model_t=None):

    # optimizer
    mpar = model.parameters()
    optimizer = torch.optim.Adam(mpar, lr=opt.learning_rate, weight_decay=0.0005)

    if torch.cuda.is_available():
        if opt.n_gpu > 1:
            model = nn.DataParallel(model)
        # model = model.cuda()
        # ll.cuda()
        
        torch.backends.cudnn.benchmark = True

    # tensorboard
    logger = tb_logger.Logger(logdir=opt.tb_dir, flush_secs=2)

    # set cosine annealing scheduler
    if opt.cosine:
        eta_min = opt.learning_rate * (opt.lr_decay_rate ** 3)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, opt.epochs, eta_min, -1)

    # train
    epoch_ = 0
    save_file = os.path.join(opt.model_folder, 'ckpt_epoch_{epoch}.pth'.format(epoch=str(epoch_).zfill(3)))
    save_checkpoint(model, optimizer, save_file, epoch_)

    ckpts = os.listdir(opt.model_folder)
    ckpts.sort()
    print(ckpts)
    if 'ckpt' in ckpts[-1] and '000' not in ckpts[-1]:
        model, optimizer, epoch_ = load_checkpoint(model, os.path.join(opt.model_folder, ckpts[-1]))
    print(epoch_)

    for epoch in range(epoch_ + 1, opt.epochs + 1):
        
        if opt.cosine:
            scheduler.step()
        else:
            adjust_learning_rate(epoch, opt, optimizer)
        
        print("==> training")
        time1 = time.time()
        train_acc, train_loss = train_epoch(epoch, train_loader, model, ll, optimizer, opt)
        time2 = time.time()
        print('epoch {}, total time {:.2f}'.format(epoch, time2 - time1))

        logger.log_value('train_acc', train_acc, epoch)
        logger.log_value('train_loss', train_loss, epoch)

        test_acc, test_acc_top5, test_loss = validate(val_loader, model, ll, opt, ll.iou)

        logger.log_value('test_acc', test_acc, epoch)
        logger.log_value('test_acc_top5', test_acc_top5, epoch)
        logger.log_value('test_loss', test_loss, epoch)

        # regular saving
        if epoch % opt.save_freq == 0:
            print('==> Saving...')
            save_file = os.path.join(opt.model_folder, 'ckpt_epoch_{epoch}.pth'.format(epoch=str(epoch).zfill(3)))
            save_checkpoint(model, optimizer, save_file, epoch, opt.n_gpu)

    # save the last model
    save_file = os.path.join(opt.model_folder, '{}_last.pth'.format(opt.model))
    save_checkpoint(model, optimizer, save_file, opt.epochs, opt.n_gpu)
