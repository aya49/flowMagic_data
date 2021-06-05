## DataSampler
## author: alice yue aya43@sfu.ca
## date created: 2021-03-21
## date modified: 2021-03-21

import torch as pt
import numpy as np

class DataSampler():
    
    def __init__(self, label, n_batch, n_cls, n_per):
        self.n_batch = n_batch
        self.n_cls = n_cls
        self.n_per = n_per

        label = np.array(label)
        self.m_ind = []
        for i in range(max(label) + 1):
            ind = np.argwhere(label == i).reshape(-1)
            ind = pt.from_numpy(ind)
            self.m_ind.append(ind)

    def __len__(self):
        return self.n_batch
    
    def __iter__(self):
        for i_batch in range(self.n_batch):
            batch = []
            classes = pt.randperm(len(self.m_ind))[:self.n_cls]
            for c in classes:
                l = self.m_ind[c]
                pos = pt.randperm(len(l))[:self.n_per]
                batch.append(l[pos])
            batch = pt.stack(batch).t().reshape(-1)
            yield batch
