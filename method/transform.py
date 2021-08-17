import torchvision.transforms as transforms
import torchvision.transforms.functional as transfunc

# for training with training samples
def transform_A(x, y, s=None):
    # Random rotation
    deg = 90 * (transforms.RandomRotation((0,359)).get_params() // 90)
    if deg != 0:
        x = transfunc.rotate(x, deg)
        y = transfunc.rotate(y, deg)

    # Random horizontal flipping
    if random.random() > 0.5:
        x = transfunc.hflip(x)
        y = transfunc.hflip(y)

    # Random vertical flipping
    if random.random() > 0.5:
        x = transfunc.vflip(x)
        y = transfunc.vflip(y)
    
    # Random resize
    params = transforms.RandomResizedCrop(size=s, scale=(0.8, 1.0), ratio=(0.75, 1.33)).get_params()
    x = transfunc.crop(x, *params)
    y = transfunc.crop(x, *params)
    
    return x, y

# for meta training with reference samples
def transform_B(x, y, s):
    # Random resize
    params = transforms.RandomResizedCrop(size=s, scale=(0.8, 1.0), ratio=(0.75, 1.33)).get_params()
    x = transfunc.crop(x, *params)
    y = transfunc.crop(x, *params)

    return x, y


transforms_dict = {
    'A': transform_A,
    'B': transform_B
}


transform_names = list()
for name, dict_ in transforms_dict.items():
    transform_names.append(name)