import numpy as np
from skimage import io, transform
from scipy.ndimage import gaussian_filter,  gaussian_gradient_magnitude
import time
import matplotlib.pyplot as plt

d=time.time()

image_complete = io.imread("resultat_retouche.png")
image_complete = image_complete[:,:,:3] #enlever la composante de transparence

def filtres(image, outputres):
    image=transform.resize(image_complete,outputres,anti_aliasing=True)
    imageLowRes=transform.resize(image_complete,(270,480),anti_aliasing=True)

    #Masques a ajouter en basse résolution
    f=1920./imageLowRes.shape[0]
    print("début filtre ",time.time()-d)
    f1=time.time()
    coeff_sat = 0.3
    saturation = np.copy(imageLowRes)
    saturation[saturation>coeff_sat]=1
    saturation[saturation<=coeff_sat]=0
    print("saturation de l'image ", time.time()-f1)
    f1=time.time()
    contours= gaussian_gradient_magnitude(saturation, 250/f)
    contours[:,:,1]=0.4*contours[:,:,1]
    contours[:,:,2]=0.
    print("gradient gaussien ", time.time()-f1)
    f1=time.time()
    contours2 = gaussian_filter(saturation, sigma = 500/f)
    contours2[:,:,1]=0.8*contours2[:,:,1]
    contours2[:,:,2]=0.3*contours2[:,:,2]
    print("gaussien ",time.time()-f1)
    f1=time.time()

    #Les mettres à la même resolution que l'image de sortie
    
    contoursup=transform.resize(contours,outputres,anti_aliasing=True)
    contours2up=transform.resize(contours2,outputres,anti_aliasing=True)

    image_filtree = image+2*contours2up+30/f*contoursup
    image_filtree[image_filtree>1]=1
    print("saturation de l'image ",time.time()-f1)
    return image_filtree

image_filtree = filtres(image_complete,(1080,1920))

io.imsave("image_filtree.png",image_filtree)

f=time.time()
print(f"Temps d'exécution de {f-d} s ")