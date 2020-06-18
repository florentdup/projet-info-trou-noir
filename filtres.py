import numpy as np
import skimage
from skimage import io, transform
import scipy.optimize as op
from scipy.ndimage import percentile_filter, gaussian_filter,  gaussian_gradient_magnitude
import time
import matplotlib.pyplot as plt

image_complete = io.imread("resultat_retouche.png")
image_complete = image_complete[:,:,0:3]

def filtres(image, reduction=5):
    image = skimage.transform.resize(image, (image.shape[0]/reduction,image.shape[1]/reduction,image.shape[2]))
    image = gaussian_filter(image, sigma=1)
    contours= gaussian_gradient_magnitude(image, 250/reduction)
    contours[:,:,1]=0.4*contours[:,:,1]
    contours[:,:,2]=0*contours[:,:,2]
    io.imsave("contours.png",contours)
    print('gradient gaussien : ok')
    pourcent = percentile_filter(image, percentile=20, size = int(100/reduction))
    io.imsave("pourcent.png",pourcent)
    print('filtre pourcent : ok')
    contours2 = gaussian_filter(pourcent, sigma = 500/reduction)
    contours2[:,:,1]=0.8*contours2[:,:,1]
    contours2[:,:,2]=0.3*contours2[:,:,2]
    io.imsave("contours2.png",contours2)
    image_filtree = image+3*contours2+300/reduction*contours
    image_filtree[image_filtree>1]=1
    return image_filtree

image_filtree = filtres(image_complete)
print(image_filtree.shape)
plt.imshow(image_filtree)

io.imsave("image_filtree_kerr_4.png",image_filtree)