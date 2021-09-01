import torch
import torch.nn as nn

from mmseg.models.segmentors.base import BaseSegmentor

from .vit_mla_1backbone import VIT_MLA_T
from .vit_mla_2head import VIT_MLAHead
from .vit_mla_3auxihead import VIT_MLA_AUXIHead

class EncoderDecoder(BaseSegmentor):
    
    def __init__(self, opt, init_path=None,
        
        # # VIT_MLA
        # img_size=200, patch_size=10, pos_embed_interp=True, drop_rate=0.,
        # depth=12, patch_size=10, in_chans=2, 
        # # self, model_name='vit_large_patch16_384', 
        # # img_size=384, patch_size=16, 
        # # in_chans=3, embed_dim=1024, depth=24,
        # # num_heads=16, # num_classes=19, 
        # # mlp_ratio=4., qkv_bias=True, qk_scale=None, 
        # # drop_rate=0.1, attn_drop_rate=0., drop_path_rate=0., 
        # # hybrid_backbone=None, 
        # # norm_layer=partial(nn.LayerNorm, eps=1e-6), norm_cfg=None,
        # # pos_embed_interp=False, random_init=False, align_corners=False, 
        # # mla_channels=256, mla_index=(5, 11, 17, 23), 
        
        **kwargs):
        
        super(EncoderDecoder, self).__init__()
        
        self.backbone = VIT_MLA_T(
            img_size=opt.dim, patch_size=20,
            in_chans=3, embed_dim=1024, depth=opt.depth,
            num_heads=16, num_classes=opt.n_class,
            mlp_ratio=4., 
            qkv_bias=True, qk_scale=None,
            pos_embed_interp=False 
            # random_init=opt.mode == 'pretrain'
        )
        
        self.decode_head = VIT_MLAHead(
            mla_channels=256, mlahead_channels=128
        )
        
        # if opt.mode == 'pretrain':
        #     self.auxiliary_head = nn.ModuleList([
        #         VIT_MLA_AUXIHead(
        #             in_channels=256,
        #             channels=512,
        #             num_classes=opt.n_class,
        #             align_corners=False
        #             in_index=i,
        #             img_size=opt.img_size,
        #             loss_decode=dict(type='LovaszLoss', loss_weight=0.4))
        #         for i in range(4)])
    
    def forward(self, x):
        
        # run the shared layer(s)
        x = self.backbone(x)

        # run the different heads with the output of the shared layers as input
        x = self.decode_head(x)
        # if opt.mode == 'pretrain':
        #     linear_out = self.auxiliary_head(x)

        return x

def model_setr(opt):
    return EncoderDecoder(opt)