import torchvision.transforms as tr
import torchvision.transforms.functional as tf

# for training with training samples
def transform_A(x, y, s=None):
    # Random rotation
    deg = 90 * (tr.RandomRotation((0,359)).get_params() // 90)
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
    params = tr.RandomResizedCrop(size=s, scale=(0.8, 1.0), ratio=(0.75, 1.33)).get_params()
    x = tf.crop(x, *params)
    y = tf.crop(x, *params)
    
    return x, y

# for meta training with reference samples
def transform_B(x, y, s):
    # Random resize
    params = tr.RandomResizedCrop(size=s, scale=(0.8, 1.0), ratio=(0.75, 1.33)).get_params()
    x = tf.crop(x, *params)
    y = tf.crop(x, *params)

    return x, y


transform_dict = {
    'A': transform_A,
    'B': transform_B
}

transform_names = list()
for name, dict_ in transform_dict.items():
    transform_names.append(name)