_base_ = '../_base_/models/setr_mla.py'
model = dict(
    type='EncoderDecoder',
    backbone=dict(
        type='VisionTransformer',
        img_size=(200, 200), # 480
        patch_size=16,
        in_channels=2, # 3
        embed_dims=256, # 512
        num_layers=6, # 24
        num_heads=16,
        out_indices=(5, 11, 17, 23),
        drop_rate=0.1,
        # norm_cfg=backbone_norm_cfg,
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
        img_size=(200, 200), #
        in_channels=(128, 128, 128, 128), # 256
        channels=256, # 512
        in_index=(0, 1, 2, 3),
        dropout_ratio=0,
        mla_channels=64, # 128
        num_classes=5, ###
        align_corners=False,
        loss_decode=dict(type='LovaszLoss', loss_type='multi_class')
    )
)