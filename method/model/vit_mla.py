import mmcv
from mmseg.apis import init_segmentor, inference_segmentor, init_cfg

cfg = mmcv.Config.fromfile('vit_mla_cfg.py')
# cfg = mmcv.Config.fromfile('/home/aya43/flowMagic_data/src/method/model/vit_mla_cfg.py')
model = init_segmentor(cfg)