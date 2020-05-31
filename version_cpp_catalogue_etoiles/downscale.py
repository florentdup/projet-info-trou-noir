from PIL import Image,ImageFilter
import matplotlib.pyplot as plt

img = Image.open("resultat.png")

width, height = img.size
factor=0.5

img_upscaled = img.resize((int(width*factor),int(height*factor)),resample=Image.ANTIALIAS)

img_upscaled.save("resultat_FHD.png")
