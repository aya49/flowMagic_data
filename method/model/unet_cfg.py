norm_cfg = dict(type='BN', requires_grad=True)
model = dict(
    type='EncoderDecoder',
    pretrained=None,
    backbone=dict(
        type='UNet',
        in_channels=2, ##
        base_channels=32,
        num_stages=4,
        strides=(1, 1, 1, 1),
        enc_num_convs=(2, 2, 2, 2),
        dec_num_convs=(2, 2, 2),
        downsamples=(True, True, True),
        enc_dilations=(1, 1, 1, 1),
        dec_dilations=(1, 1, 1),
        with_cp=False,
        conv_cfg=None,
        norm_cfg=norm_cfg,
        act_cfg=dict(type='ReLU'),
        upsample_cfg=dict(type='InterpConv'),
        norm_eval=False
    ),
    decode_head=dict(
        type='ASPPHead',
        in_channels=32,
        in_index=4,
        channels=8,
        dilations=(1, 6, 12, 16),
        dropout_ratio=0.1,
        num_classes=6, ##
        norm_cfg=norm_cfg,
        align_corners=False,
        loss_decode=dict(type='LovaszLoss', loss_type='multi_class', per_image=True)
    ),
    train_cfg=dict(mode='whole'),
    test_cfg=dict(mode='whole')
)
evaluation = dict(metric='mIoU')