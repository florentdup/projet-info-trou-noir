import kerr


## Interface utilisateur trou noirs

from tkinter import *
from tkinter import ttk
from tkinter import Toplevel
from PIL import ImageTk, Image
import numpy as np


import numpy as np
from skimage import io, transform
from scipy.ndimage import gaussian_filter,  gaussian_gradient_magnitude
import matplotlib.pyplot as plt

def filtres(image, outputres):
    image=transform.resize(image,outputres,anti_aliasing=True)
    imageLowRes=transform.resize(image,(270,480),anti_aliasing=True)



    #Masques a ajouter en basse résolution
    f=1920./imageLowRes.shape[0]
    
    coeff_sat = 0.3
    saturation = np.copy(imageLowRes)
    saturation[saturation>coeff_sat]=1
    saturation[saturation<=coeff_sat]=0
    print("saturation de l'image")
    
    contours= gaussian_gradient_magnitude(saturation, 250/f)
    contours[:,:,1]=0.4*contours[:,:,1]
    contours[:,:,2]=0.
    print("gradient gaussien")
    
    contours2 = gaussian_filter(saturation, sigma = 500/f)
    contours2[:,:,1]=0.8*contours2[:,:,1]
    contours2[:,:,2]=0.3*contours2[:,:,2]
    print("gaussien")
    

    #Les mettres à la même resolution que l'image de sortie
    
    contoursup=transform.resize(contours,outputres,anti_aliasing=True)
    contours2up=transform.resize(contours2,outputres,anti_aliasing=True)


    

    image_filtree = image+2*contours2up+150/f*contoursup
    image_filtree[image_filtree>1]=1
    print("saturation de l'image")
    return image_filtree

## page d'accueil

fenetre = Tk()
fenetre.title('projet trou noir')
fenetre.configure(background="black")

label1 = Label(fenetre, text="Bienvenue dans le simulateur \n de trous noirs !", font=("Franklin gothic heavy", 22), bg='black', fg='white')
label1.pack(pady=10,padx=20)

label2 = Label(fenetre, text="par Florent Dupont, Marie-Clémentine Quilleriet \n et Camille Srecki", font=("Calibri", 11), bg='black', fg='white')
label2.pack(padx=10, pady=5)

surech=False

def start_window():
    params = Toplevel(fenetre,padx=5,pady=5)
    params.title = ('choix des parametres')

    label = Label(params, text='Choisissez les paramètres de la simulation :', font=("Calibri Bold", 12))
    label.grid(row=0,column=0,stick=E)

    label_gauche = Label(params, anchor=W, text='Paramètres du trou noir :', font=("Calibri", 12))
    label_gauche.grid(row=1,column=0,sticky=W)

    label_droit = Label(params, anchor=E, text='Paramètres de la caméra :', font=("Calibri", 12))
    label_droit.grid(row=1,column=1,sticky=W)

    # Parametre de Kerr
    var_a = DoubleVar()
    a = Scale(params,from_=0.05,to=0.5,digits=2, variable = var_a, resolution=0.01, orient=HORIZONTAL, label="Parametre de Kerr",length=200)
    a.set(.5)
    a.grid(row=2,rowspan=2,column=0,sticky=W)

    # Rayon maximal disque d'accrétion
    RAdisk = Scale(params,from_=3,to=25,orient=HORIZONTAL,label="Rayon maximal du disque d'accrétion",length=200,resolution=0.01)
    RAdisk.grid(row=4,column=0,sticky=W)
    RAdisk.set(16.)

    # position de la caméra
    camera_pos_x = Scale(params,from_=-25,to=25,orient=HORIZONTAL,label="Coordonnée x de la caméra",length=200,resolution=0.01)
    camera_pos_x.grid(row=4,column=1,sticky=W)
    camera_pos_x.set(20)


    camera_pos_y = Scale(params,from_=-25,to=25,orient=HORIZONTAL,label="Coordonnée y de la caméra",length=200,resolution=0.01)
    camera_pos_y.grid(row=5,column=1,sticky=W)
    camera_pos_y.set(0)

    camera_pos_z = Scale(params,from_=-25,to=25,orient=HORIZONTAL,label="Coordonnée z de la caméra",length=200,resolution=0.01)
    camera_pos_z.grid(row=6,column=1,sticky=W)
    camera_pos_z.set(1.2)

    # Rayon à l'infini : il faut implémenter un callback ou une binding methode

    d = [camera_pos_x.get(), camera_pos_y.get(), camera_pos_z.get()]
    d = np.array(d)
    d = np.linalg.norm(d)

    R_inf = Scale(params,from_= max(15.,RAdisk.get(),d), to=40,orient=HORIZONTAL,label="Rayon à l'infini",length=200,resolution=0.01)
    R_inf.grid(row=5,column=0,sticky=W)

    checkbuttonsurech=ttk.Checkbutton(params, text="sur-échantillonage")
    checkbuttonsurech.grid(row=6,column=0, sticky=W)
    checkbuttonsurech.state(['!alternate'],)
    checkbuttonsurech.state(['!selected'])


    # résolution
    listres = [f"240,135",f"480,270",f"1920,1080"]
    
    variable = StringVar(params)
    variable.set(listres[2])
    res = OptionMenu(params,variable,*listres)
    label = Label(params,text="Résolution du rendu")
    res.grid(row=3,column=1)
    label.grid(row=2,column=1)


    #bouton de lancement
    def demarrer_programme():

        surech=checkbuttonsurech.instate(['selected'])


        tmp=variable.get()[:].split(',') #Faire une copie de la string, puis la séparer 

        paramresolution="" 

        height=int(tmp[0])
        width=int(tmp[1])
        
        if surech: 
            paramresolution=str(2*height)+' '+str(2*width) 
        else:
            paramresolution=str(height)+' '+str(width) 

    
        liste_param = open("params.txt",'w')
        liste_param.write(f"{a.get()} {RAdisk.get()} {camera_pos_x.get()} {camera_pos_y.get()} {camera_pos_z.get()} {R_inf.get()} {paramresolution}")
        liste_param.close()
        kerr.image("params.txt")

        #Charger le resultat
        image = io.imread("resultat.png")

        imageFinale=transform.resize(image,(width,height),anti_aliasing=True)

        io.imsave("resultatFinal.png",imageFinale)

        imagefinaleFiltree = filtres(imageFinale,(width,height))
        io.imsave("resultatFinalFiltré.png",imagefinaleFiltree)

        io.imshow(imageFinale)
        plt.show()
        io.imshow(imagefinaleFiltree)
        plt.show()



        print("ok")
        

    bouton_lancer = Button(params, text = 'Démarrer', command=demarrer_programme, fg='red', activebackground='white')
    bouton_lancer.grid(row=8, columnspan = 2, pady=10 )

    params.mainloop()


bouton_commencer = Button(fenetre, text = 'Commencer', command=start_window, bg='black', fg='orange', activebackground='white')
bouton_commencer.pack(pady=5 )

img = ImageTk.PhotoImage(Image.open("image_filtree_kerr_3.png"))
panel = Label(fenetre, image = img,bg='black')
panel.pack(ipady=0,ipadx=0,padx=20)

label3 = Label(fenetre, text="Projet réalisé dans le cadre de l'UE d'informatique \n du cycle ingénieur civil de Mines ParisTech \n sous la tutelle de Nikolas Stott", font=("Calibri", 8), bg='black', fg='white')
label3.pack(padx=10,pady=10)

#bouton_biblio = Button(fenetre, text = 'Sources', command=start_window, bg='black', fg='white', activebackground='white')
#bouton_biblio.pack(pady=10)

fenetre.mainloop()



