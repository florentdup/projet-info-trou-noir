from PIL import Image,ImageFilter
import matplotlib.pyplot as plt

img = Image.open("adisk.png")

width, height = img.size
factor=4
img=img.filter(ImageFilter.GaussianBlur(3))


#plt.imshow(img)
#plt.show()


img_upscaled = img.resize((width*factor,height*factor),resample=Image.BILINEAR)

#plt.imshow(img_upscaled)
#plt.show()

img_upscaled.save("adisk_upscaled.png")

