_base_ = '../_base_/models/setr_mla.py'
model = dict(
    backbone=dict(img_size=opt.dim, pos_embed_interp=True, drop_rate=0.,
                  mla_channels=128, mla_index=(5, 8, 17, 23)),
    decode_head=dict(img_size=opt.dim, mla_channels=128,
                     mlahead_channels=64, num_classes=5)
)