import torch
import mmcv
from mmseg.apis import init_segmentor#, inference_segmentor, init_cfg

def model_unet(opts):
    cfg = mmcv.Config.fromfile('unet_cfg.py')
    # cfg = mmcv.Config.fromfile('/home/aya43/flowMagic_data/src/method/model/vit_mla_cfg.py')
    
    if torch.cuda.is_available():
        model = init_segmentor(cfg, device='cuda:0')
    else:
        model = init_segmentor(cfg, device='cpu')

    return model