import torchvision.transforms as tr
import torchvision.transforms.functional as tf

import random

tresize = tr.Compose([tr.Resize((256, 256))])
 
# for training with training samples
def transform_A(x, y):
    # Random rotation
    deg = tr.RandomRotation.get_params(degrees=(0,45))
    if random.random() > 0.5:
        deg = 360-deg

    if deg != 0:
        x = tf.rotate(x, deg)
        y = tf.rotate(y, deg)

    # Random horizontal flipping
    if random.random() > 0.5:
        x = tf.hflip(x)
        y = tf.hflip(y)
    
    # Random vertical flipping
    if random.random() > 0.5:
        x = tf.vflip(x)
        y = tf.vflip(y)
    
    # Random resize
    params = tr.RandomResizedCrop.get_params(img=x, scale=(0.8, 1.0), ratio=(0.75, 1.33))
    x = tf.crop(x, *params)
    y = tf.crop(y, *params)

    x = tresize(x)
    y = tresize(y)
    
    return x, y

# for meta training with reference samples
def transform_B(x, y):
    # Random resize
    params = tr.RandomResizedCrop.get_params(img=x, scale=(0.8, 1.0), ratio=(0.75, 1.33))
    x = tf.crop(x, *params)
    y = tf.crop(y, *params)

    x = tresize(x)
    y = tresize(y)

    return x, y


transform_dict = {
    'A': transform_A,
    'B': transform_B
}

transform_names = list()
for name, dict_ in transform_dict.items():
    transform_names.append(name)