## Interface utilisateur trou noirs

from tkinter import *
from tkinter import ttk
from tkinter import Toplevel
from PIL import ImageTk, Image
import numpy as np

## page d'accueil

fenetre = Tk()
fenetre.title('projet trou noir')
fenetre.configure(background="black")

label1 = Label(fenetre, text="Bienvenue dans le simulateur \n de trous noirs !", font=("Franklin gothic heavy", 22), bg='black', fg='white')
label1.pack(pady=10,padx=20)

label2 = Label(fenetre, text="par Florent Dupont, Marie-Clémentine Quilleriet \n et Camille Srecki", font=("Calibri", 11), bg='black', fg='white')
label2.pack(padx=10, pady=5)

def start_window():
    params = Toplevel(padx=5,pady=5)
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
    a.grid(row=2,rowspan=2,column=0,sticky=W)

    # Rayon du disque d'accrétion
    RAdisk = Scale(params,from_=2,to=6,orient=HORIZONTAL,label="Rayon du disque d'accrétion",length=200)
    RAdisk.grid(row=4,column=0,sticky=W)

    # position de la caméra
    camera_pos_x = Scale(params,from_=-20,to=20,orient=HORIZONTAL,label="Coordonnée x de la caméra",length=200)
    camera_pos_x.grid(row=4,column=1,sticky=W)

    camera_pos_y = Scale(params,from_=-20,to=20,orient=HORIZONTAL,label="Coordonnée y de la caméra",length=200)
    camera_pos_y.grid(row=5,column=1,sticky=W)

    camera_pos_z = Scale(params,from_=-20,to=20,orient=HORIZONTAL,label="Coordonnée z de la caméra",length=200)
    camera_pos_z.grid(row=6,column=1,sticky=W)

    # Rayon max du disque
    R_max = Scale(params,from_=8,to=32,orient=HORIZONTAL,label="Rayon maximal du disque",length=200)
    R_max.grid(row=5,column=0,sticky=W)
    print(type(R_max))

    d = [camera_pos_x.get(), camera_pos_y.get(), camera_pos_z.get()]
    d = np.array(d)
    d = np.linalg.norm(d)

    # Rayon à l'infini : il faut implémenter un callback ou une binding methode
    R_inf = Scale(params,from_= max(R_max.get(),d), to=100,orient=HORIZONTAL,label="Rayon à l'infini",length=200)
    R_inf.grid(row=6,column=0,sticky=W)

    # résolution
    listres = [f"384,216",f"960,540",f"1920,1080"]
    variable = StringVar(params)
    variable.set(listres[2])
    res = OptionMenu(params,variable,*listres)
    label = Label(params,text="Résolution du rendu")
    res.grid(row=3,column=1)
    label.grid(row=2,column=1)

    # Mode : grille, photo, photo filtree


    #tkMessageBox si erreur
    #bouton de lancement def ...


bouton_commencer = Button(fenetre, text = 'Commencer', command=start_window, bg='black', fg='orange', activebackground='white')
bouton_commencer.pack(pady=5 )


img = ImageTk.PhotoImage(Image.open("C:\\Users\\camil\\Desktop\\Mines 2019-2020\\Info\\projet\\image_filtree_kerr_3.png"))
panel = Label(fenetre, image = img,bg='black')
panel.pack(ipady=0,ipadx=0,padx=20)

label3 = Label(fenetre, text="Projet réalisé dans le cadre de l'UE d'informatique \n du cycle ingénieur civil de Mines ParisTech \n sous la tutelle de Nikolas Stott", font=("Calibri", 8), bg='black', fg='white')
label3.pack(padx=10,pady=10)

bouton_biblio = Button(fenetre, text = 'Sources', command=start_window, bg='black', fg='white', activebackground='white')
bouton_biblio.pack(pady=10)

fenetre.mainloop()