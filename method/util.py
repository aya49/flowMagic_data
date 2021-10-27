import torch
import time
import numpy as np

from matplotlib import pyplot as plt

def visualize(save_file, **images):
    n = len(images)
    plt.figure(figsize=(16, 5))
    for i, (name, image) in enumerate(images.items()):
        plt.subplot(1, n, i + 1)
        plt.xticks([])
        plt.yticks([])
        plt.title(' '.join(name.split('_')).title())
        plt.imshow(image)
    plt.savefig(save_file)

def prep_input(inp, target, xfn):
    if torch.cuda.is_available():
        inp = inp.cuda()
        target = target.cuda()
    
    (H, W, C) = (inp.shape[3], inp.shape[2], inp.shape[1])
    img_metas = [{
        'img_shape': (H, W, C),
        'ori_shape': (H, W, C),
        'pad_shape': (H, W, C),
        'filename': xfn_,
        'scale_factor': 1.0,
        'flip': False,
    } for xfn_ in xfn]

    return inp, target, img_metas

def save_checkpoint(model, optimizer, save_path, epoch, n_gpu=1):
    torch.save({
        'model': model.state_dict() if n_gpu <= 1 else model.module.state_dict(),
        'optimizer': optimizer.state_dict(),
        'epoch': epoch
    }, save_path)

def load_checkpoint(model, save_path):
    model_state = torch.load(save_path)
    model.load_state_dict(model_state['model'])
    return model, model_state['optimizer'], model_state['epoch']

def get_teacher_name(model_dir):
    """parse to get teacher model name"""
    segments = model_dir.split('/')[-2].split('_')
    if ':' in segments[0]:
        return segments[0].split(':')[-1]
    else:
        if segments[0] != 'wrn':
            return segments[0]
        else:
            return segments[0] + '_' + segments[1] + '_' + segments[2]

def adjust_learning_rate(epoch, opt, optimizer):
    """Sets the learning rate to the initial LR decayed by decay rate every steep step"""
    steps = np.sum(epoch > np.asarray(opt.lr_decay_epochs))
    if steps > 0:
        new_lr = opt.learning_rate * (opt.lr_decay_rate ** steps)
        for param_group in optimizer.param_groups:
            param_group['lr'] = new_lr
            
class AverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0
        self.vals = []

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count
        self.vals.append(val)
