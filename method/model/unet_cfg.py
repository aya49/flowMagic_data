norm_cfg = dict(type='BN', requires_grad=True)
model = dict(
    type='EncoderDecoder',
    pretrained=None,
    backbone=dict(
        type='UNet',
        in_channels=2, ##
        base_channels=16,
        num_stages=5,
        strides=(1, 1, 1, 1, 1),
        enc_num_convs=(2, 2, 2, 2, 2),
        dec_num_convs=(2, 2, 2, 2),
        downsamples=(True, True, True, True),
        enc_dilations=(1, 1, 1, 1, 1),
        dec_dilations=(1, 1, 1, 1),
        with_cp=False,
        conv_cfg=None,
        norm_cfg=norm_cfg,
        act_cfg=dict(type='ReLU'),
        upsample_cfg=dict(type='InterpConv'),
        norm_eval=False
    ),
    decode_head=dict(
        type='ASPPHead',
        in_channels=16,
        in_index=4,
        channels=8,
        dilations=(1, 12, 24, 36),
        dropout_ratio=0.1,
        num_classes=5, ##
        norm_cfg=norm_cfg,
        align_corners=False,
        loss_decode=dict(type='LovaszLoss', loss_type='multi_class', per_image=True)
    ),
    train_cfg=dict(crop_size=(256, 256)), ##
    test_cfg=dict(mode='slide', crop_size=(256, 256), stride=(170, 170)) ##
)
evaluation = dict(metric='mIoU')