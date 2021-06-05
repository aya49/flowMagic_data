# generator functions that load 2D scatterplot matrices and ground truth T/F matrices
# written 2020-11
# updated 2020-12
# by alice yue

# imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import sklearn.datasets as sd
import gen_func as gf
import matplotlib.pyplot as plt

import os
import glob
from pathlib import Path

import re

from sklearn.datasets.samples_generator import make_blobs

golden_path = "/mnt/f/Brinkman group/COVID/data/structure_test/golden_samples"
data_path = "/mnt/f/Brinkman group/COVID/data/structure_test/data"

golden_files = np.array(Path(golden_path).rglob("*.csv"))
data_files = [Path(re.sub(golden_path, data_path, str(x))) for x in golden_files]
data_exists = [os.path.exists(x) for x in data_files]

n_components = 3
X, truth = sd.make_blobs(n_samples=300, centers=n_components,
                      cluster_std = [2, 1.5, 1], random_state=42)

kl, kg = gf.kde2D(X, xbins=100j, ybins=100j)
