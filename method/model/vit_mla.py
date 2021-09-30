import mmcv
from mmseg.apis import init_segmentor, inference_segmentor

cfg = mmcv.Config.fromfile('vit_mla_cfg.py')