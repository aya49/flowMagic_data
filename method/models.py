from genericpath import samefile
import torch
import mmcv
from mmseg.apis import init_segmentor#, inference_segmentor, init_cfg
from mmseg.models import build_segmentor
from mmcv import ConfigDict

import segmentation_models_pytorch as smp

def model_unet(opt):
    model = smp.Unet(
        encoder_name="resnet18",        # encoder
        encoder_depth=5,
        # encoder_weights="None",       # random initialization
        in_channels=len(opt.x_2D),      # model input channels (1 for gray-scale images, 3 for RGB, etc.)
        classes=opt.n_class,            # model output channels (number of classes in your dataset)
    )
    return model

def model_unet_(opt):
    cfg = mmcv.Config.fromfile('model/unet_cfg.py')
    # cfg = mmcv.Config.fromfile('/home/aya43/flowMagic_data/src/method/model/unet_cfg.py')
    
    model = init_segmentor(cfg)
    return model

def model_setr(opt):
    cfg = mmcv.Config.fromfile('model/vit_mla_cfg.py')
    # cfg = mmcv.Config.fromfile('/home/aya43/flowMagic_data/src/method/model/vit_mla_cfg.py')
    
    model = init_segmentor(cfg)
    return model

model_dict = {
    'setr': model_setr,
    'unet': model_unet
}

model_names = list()
for name, dict_ in model_dict.items():
    model_names.append(name)

def create_model(opt):
    return model_dict[opt.model](opt)



def metafreeze_model(model, opt):
    # freeze
    for param in model.parameters():
        param.requires_grad = False
    
    if opt.model == 'setr':
        for param in model.backbone.layers[6].parameters():
            param.requires_grad = True
        for param in model.backbone.layers[7].parameters():
            param.requires_grad = True
    
    if opt.model == 'unet':
        for param in model.backbone.encoder[3].parameters():
            param.requires_grad = True
        for param in model.backbone.decoder[0].parameters():
            param.requires_grad = True
    
    return model
        