import torch
import mmcv
from mmseg.apis import init_segmentor#, inference_segmentor, init_cfg
from mmcv import ConfigDict

import segmentation_models_pytorch as smp

def model_unet_(opt):
    model = smp.Unet(
        encoder_name="resnet18",        # encoder
        encoder_depth=opt.depth,
        # encoder_weights="None",       # random initialization
        in_channels=len(opt.x_2D),      # model input channels (1 for gray-scale images, 3 for RGB, etc.)
        classes=opt.n_class,            # model output channels (number of classes in your dataset)
    )

    return model

def model_unet(opt):
    # cfg = mmcv.Config.fromfile('/home/aya43/flowMagic_data/src/method/model/unet_cfg.py')
    
    cfg = ConfigDict(
        type='EncoderDecoder',
        pretrained=None,
        backbone=dict(
            type='UNet',
            in_channels=len(opt.x_2D),
            base_channels=64,
            num_stages=5,
            strides=(1, 1, 1, 1, 1),
            enc_num_convs=(2, 2, 2, 2, 2),
            dec_num_convs=(2, 2, 2, 2),
            downsamples=(True, True, True, True),
            enc_dilations=(1, 1, 1, 1, 1),
            dec_dilations=(1, 1, 1, 1),
            with_cp=False,
            conv_cfg=None,
            norm_cfg=dict(type='BN', requires_grad=True),
            act_cfg=dict(type='ReLU'),
            upsample_cfg=dict(type='InterpConv'),
            norm_eval=False),
        decode_head=dict(
            type='ASPPHead',
            in_channels=64,
            in_index=4,
            channels=16,
            dilations=(1, 12, 24, 36),
            dropout_ratio=0.1,
            num_classes=2,
            norm_cfg=dict(type='BN', requires_grad=True),
            align_corners=False,
            loss_decode=dict(type='LovaszLoss', loss_type='multi_class', per_image=True)
        ),
        # model training and testing settings
        train_cfg=dict(crop_size=opt.dim),
        test_cfg=dict(mode='slide', crop_size=opt.dim, stride=170)
    )

    if torch.cuda.is_available():
        model = init_segmentor(cfg, device='cuda:0')
    else:
        model = init_segmentor(cfg, device='cpu')

    return model

def model_setr(opt):
    # cfg = mmcv.Config.fromfile('/home/aya43/flowMagic_data/src/method/model/vit_mla_cfg.py')
    cfg =  = ConfigDict(
    type='EncoderDecoder',
        backbone=dict(
            type='VisionTransformer',
            img_size=(opt.dim, opt.dim), # 480
            patch_size=16,
            in_channels=2, # 3
            embed_dims=opt.dim, # 512
            num_layers=6, # 24
            num_heads=16,
            out_indices=(5, 11, 17, 23),
            drop_rate=0.1,
            norm_cfg=dict(type='BN', requires_grad=True),
            with_cls_token=False,
            interpolate_mode='bilinear',
        ),
        neck=dict(
            type='MLANeck',
            in_channels=[512, 512, 512, 512], # 1024
            out_channels=128, # 256
            act_cfg=dict(type='ReLU'),
        ),
        decode_head=dict(
            type='SETRMLAHead',
            in_channels=(128, 128, 128, 128), # 256
            channels=opt.dim, # 512
            in_index=(0, 1, 2, 3),
            dropout_ratio=0,
            mla_channels=64, # 128
            num_classes=5, ###
            align_corners=False,
            loss_decode=dict(type='LovaszLoss', loss_type='multi_class', per_image=True)
        ),
        train_cfg=dict(crop_size=opt.dim),
        test_cfg=dict(mode='slide', crop_size=opt.dim, stride=170)
    )
    
    if torch.cuda.is_available():
        model = init_segmentor(cfg, device='cuda:0')
    else:
        model = init_segmentor(cfg, device='cpu')

    return model

model_dict = {
    'setr': model_setr,
    'unet': model_unet
}

model_names = list()
for name, dict_ in model_dict.items():
    model_names.append(name)

def create_model(opt):
    return model_dict[name](opt)
        