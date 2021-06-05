import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib.colors as plc
import matplotlib.image as img
import mpl_scatter_density # adds projection='scatter_density'
import cv2

# outputs the (100x100)x3 point cloud and 100x100 kernel density grid given nx2 matrix
def kde2D(X, Y, xbins=100j, ybins=100j, png_name="temp.png", dpi=100, **kwargs):
    x = X[:, 0]
    y = X[:, 1]

    # create grid 100x100
    xx, yy = np.mgrid[x.min():x.max():xbins, y.min():y.max():ybins]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])

    # 2D kernel density
    kernel = st.gaussian_kde(values, **kwargs)

    # out: grid
    kg = np.reshape(kernel(positions).T, xx.shape)

    # out: point cloud
    kgf = kg.flatten()
    kl = np.zeros((len(kgf), 3))
    kl[:, 0] = xx.flatten()
    kl[:, 1] = yy.flatten()
    kl[:, 2] = kg.flatten()/max(kg)

    if (png_name != None):
        greyscale = plc.LinearSegmentedColormap.from_list(
            'greyscale', [(0, '#ffffff'), (1, '#000000')], N=256)
        fig = plt.figure(figsize=(xx.shape[0]//dpi, xx.shape[1]//dpi), dpi=dpi)
        ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
        density = ax.scatter_density(x, y, cmap=greyscale)

        # remove axis
        plt.axis('off')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        plt.savefig(png_name, bbox_inches='tight', pad_inches=0)

        figi = img.imread(png_name)
        figg = cv2.cvtColor(figi, cv2.COLOR_RGB2GRAY)
        figm = np.array(figg)
        figm = figm/max(figm)

    return kl, kg, figm
