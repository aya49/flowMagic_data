# module load LIB/CUDA/10.2

# !pip install segmentation_models_pytorch
# !pip install torchviz
# !MMCV_WITH_OPS=1 pip install mmcv-full==1.2.2 -f https://download.openmmlab.com/mmcv/dist/cu110/torch1.7.0/index.html #full version not for windows

import sys
import os
import trace

import pandas as pd
import numpy as np

import torch
from torch.utils.data import Dataset
import torchvision.transforms as transforms
import torchvision.transforms.functional as FT
import torch.nn as nn
import torch.nn.functional as FN

from torchviz import make_dot

import argparse
import time

import mmcv
import mmseg

from mmcv.utils import Config
from mmseg.models import build_segmentor

print("cuda available")
print(torch.cuda.is_available())

cfg = Config.fromfile("configs/SETR/SETR_MLA_480x480_80k_pascal_context_bs_16.py")
model = build_segmentor(cfg.model, train_cfg=cfg.train_cfg, test_cfg=cfg.test_cfg)

