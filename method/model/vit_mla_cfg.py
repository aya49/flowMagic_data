backbone_norm_cfg = dict(type='LN', eps=1e-6, requires_grad=True)
norm_cfg = dict(type='BN', requires_grad=True)
model = dict(
    type='EncoderDecoder',
    backbone=dict(
        type='VisionTransformer',
        img_size=(256, 256), # 480
        patch_size=8,
        in_channels=2, # 3
        embed_dims=128, # 512
        num_layers=12, # 24
        num_heads=8,
        out_indices=(3, 4, 5, 6),
        drop_rate=0.1,
        norm_cfg=backbone_norm_cfg,
        with_cls_token=False,
        interpolate_mode='bilinear',
    ),
    neck=dict(
        type='MLANeck',
        in_channels=[256, 256, 256, 256], # 1024
        out_channels=64, # 256
        act_cfg=dict(type='ReLU'),
    ),
    decode_head=dict(
        type='SETRMLAHead',
        in_channels=(64, 64, 64, 64), # 256
        channels=128, # 512
        in_index=(0, 1, 2, 3),
        dropout_ratio=0,
        mla_channels=64, # 128
        num_classes=5, ###
        align_corners=False,
        loss_decode=dict(type='LovaszLoss', loss_type='multi_class', per_image=True)
    ),
    train_cfg=dict(crop_size=256),
    test_cfg=dict(mode='slide', crop_size=256, stride=170)
)