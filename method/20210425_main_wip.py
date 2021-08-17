import torch

import mmcv
import mmseg

from mmcv.utils import Config
from mmseg.models import build_segmentor

print("cuda available")
print(torch.cuda.is_available())

cfg = Config.fromfile("configs/SETR/SETR_MLA_480x480_80k_pascal_context_bs_16.py")
model = build_segmentor(cfg.model, train_cfg=cfg.train_cfg, test_cfg=cfg.test_cfg)
