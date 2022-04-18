import torchvision.transforms as tr
import torchvision.transforms.functional as tf

import random
# import cv2
import numpy as np

from util import bound_rect

tresize = tr.Compose([tr.Resize((256, 256))])

def crop_y(x, y, cpop, dim=None, ratio=1/.85):
    if cpop==-1:
        cpop = random.randrange(1,int(y.max()))
    xind, yind, w, h = dim if not dim==None else bound_rect((y[0] == cpop).nonzero())
    if xind==None:
        return None
    w_ = int(w*ratio) # index > ncol is ok
    h_ = int(h*ratio)
    xind_ = int(max(xind-((ratio-1)/2), 0))
    yind_ = int(max(yind-((ratio-1)/2), 0))
    x = x[:,xind_:(xind_+w_),yind_:(yind_+h_)]
    y = y[:,xind_:(xind_+w_),yind_:(yind_+h_)]
    y[y!=cpop] = 0
    y[y==cpop] = 1
    
    return x, y

# for training with training samples
def transform_A(x, y, y_=None, cpop=0, dim=None, rot=True):
    if cpop>0:
        x, y = crop_y(x, y, cpop, dim)
        if x==None:
            return None
    
    if rot:
        # Random rotation
        deg = tr.RandomRotation.get_params(degrees=(0,20))
        if random.random() > 0.5:
            deg = 360-deg
        
        if deg != 0:
            x = tf.rotate(x, deg)
            y = tf.rotate(y, deg)
    
    # Random resize
    params = tr.RandomResizedCrop.get_params(img=x, scale=(0.85, 1.0), ratio=(0.85, 1.15))
    x = tf.crop(x, *params)
    y = tf.crop(y, *params)
    
    x = tresize(x)
    y = tresize(y)
        
    return x, y
    

# for test
def transform_B(x, y, y_=None, cpop=0, dim=None, rot=True):
    if cpop>0:
        x, y = crop_y(x, y, cpop, dim)
        if x==None:
            return None
    
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