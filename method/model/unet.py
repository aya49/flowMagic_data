import segmentation_models_pytorch as smp

def model_unet(opt):
    model = smp.Unet(
        encoder_name="resnet18",        # encoder
        encoder_depth=6,
        encoder_weights="None",         # random initialization
        in_channels=2,                  # model input channels (1 for gray-scale images, 3 for RGB, etc.)
        classes=opt.n_class,            # model output channels (number of classes in your dataset)
    )

    return model