## Interface utilisateur trou noirs

from tkinter import *
from tkinter import ttk
from tkinter import Toplevel
from PIL import ImageTk, Image

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

    a = Scale(params,from_=0,to=10,orient=HORIZONTAL,label="Vitesse de rotation",length=200)
    a.grid(row=2,rowspan=2,column=0,sticky=W)

    RAdisk = Scale(params,from_=2,to=6,orient=HORIZONTAL,label="Rayon du disque d'accrétion",length=200)
    RAdisk.grid(row=4,column=0,sticky=W)

    # résolution
    listres = [f"(192,108)",f"(384,216)",f"(960,540)",f"(1920,1080)"]
    variable = StringVar(params)
    variable.set(listres[2])
    res = OptionMenu(params,variable,*listres)
    label = Label(params,text="Résolution du rendu")
    res.grid(row=3,column=1)
    label.grid(row=2,column=1)


    #tkMessageBox si erreur
    #
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