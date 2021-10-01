import torch
import mmcv
from mmseg.apis import init_segmentor#, inference_segmentor, init_cfg

def model_unet(opts):
    cfg = mmcv.Config.fromfile('unet_cfg.py')
    # cfg = mmcv.Config.fromfile('/home/aya43/flowMagic_data/src/method/model/unet_cfg.py')
    
    if torch.cuda.is_available():
        model = init_segmentor(cfg, device='cuda:0')
    else:
        model = init_segmentor(cfg, device='cpu')

    return model

def model_setr(opts):
    cfg = mmcv.Config.fromfile('vit_mla_cfg.py')
    # cfg = mmcv.Config.fromfile('/home/aya43/flowMagic_data/src/method/model/vit_mla_cfg.py')
    
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
        