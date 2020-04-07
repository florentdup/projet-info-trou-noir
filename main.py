import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import time
from scipy.integrate import ode


from PIL import Image
Image.MAX_IMAGE_PIXELS = 1000000000  

start_time = time.time()
background=mpimg.imread("background.jpg")/255./1.2#eso0932a.tif
adisk=mpimg.imread("adisk.jpg")/255.
#adisk=np.flip(adisk,(0,1))#exterieur du disque en bas, sans changer les couleurs (axis 2)
#plt.imshow(adisk)
#plt.show()

print("Temps de chargement image : %s secondes " % (time.time() - start_time))

rendu="3D"

BACKGROUND_X=background.shape[0]
BACKGROUND_Y=background.shape[1]

ADISK_X=adisk.shape[0]
ADISK_Y=adisk.shape[1]
#alpha_adisk=0.95 #opacité
adisk_texture_rep=4.

STEP=1.
k1=3.**2/(20**2-1.**2)
k2=k1*(-1.)+0.4
def fstep(r2):
    return k1*r2+k2


HEIGHT=360*2 #x
WIDTH=640*2#y


#zoom sur le disque
#FOV=10
#camera_pos=np.array([-10.,1.,0.3])
#camera_dir=np.array([0.0,0.2])#theta et phi

#Vue globale
FOV=50
camera_pos=np.array([-10.,0.,0.6])
camera_dir=np.array([0.04,0.])#theta et phi


RAdisk_min=1.5
RAdisk_max=4.

RAdisk_min2=RAdisk_min**2
RAdisk_max2=RAdisk_max**2


R_schwarzschild=1. #inutile
R_min=1.*R_schwarzschild # en dessous les photons sont piégés
R_inf=20 
R_inf2=R_inf**2
R_min2=R_min**2

fov_rad=FOV*np.pi/180
distance_viewplane=WIDTH/(2*np.tan(fov_rad/2))


def RK4f(y,h2):
    f = np.zeros((6))
    f[0:3] = y[3:6]
    f[3:6] = - 1.5 * h2 * y[0:3] / sixth(y[0:3])#np.power(sqrnorm(y[0:3]),2.5) 
    return f


#def RK4f(y,L2,r6):
#    f = np.zeros(y.shape)
#    f[0:3] = y[3:6]
#    f[3:6] = - 1.5 * L2 * y[0:3] / r6
 #   return f


def normalize(x):
    return x/np.linalg.norm(x)

def sqrnorm(x):
    return np.dot(x,x)



def sixth(x):#calculer la puissance 6 de r
    tmp = sqrnorm(x)
    return tmp*tmp*tmp

def asSpherical(xyz):
    r=np.linalg.norm(xyz)
    theta   =  np.arccos(xyz[2]/r)
    phi     =  np.arctan2(xyz[1]/r,xyz[0]/r)
    return np.array([r,theta,phi])

def asCartesian(rthetaphi):
    r       = rthetaphi[0]
    theta   = rthetaphi[1]
    phi     = rthetaphi[2]
    x = r * np.sin( theta ) * np.cos( phi )
    y = r * np.sin( theta ) * np.sin( phi )
    z = r * np.cos( theta )
    return np.array([x,y,z])

#pour les vecteurs deja normés
def asSpherical_u(xyz):  
    return np.array([np.arccos(xyz[2]),np.arctan2(xyz[1],xyz[0])])

def asCartesian_u(thetaphi):
    tmp=np.sin(thetaphi[0])
    return np.array([tmp*np.cos(thetaphi[1]),tmp*np.sin(thetaphi[1]),np.cos(thetaphi[0])])



#def eq(x,L2):
#        pos=x[:3]
#        xp=x[3:]
#        F=-3./2.*L2/sixth(pos)*pos
#        return np.concatenate([xp,F])

def sim_1D(point0,angle0,maxiter=100):
    N=0
    X,Y=[],[]
    #r=np.linalg.norm(pos)
    r2=sqrnorm(point0)
    r6=0.

    y=np.zeros((6))
    
    y[0:3]=np.copy(point0)
    y[3:6]=asCartesian_u(angle0) #"vitesse"
    

    #vect=np.concatenate([point,velocity])
    L2=sqrnorm(np.cross(y[0:3],y[3:6])) #"moment cinétique"

    

    rkstep=STEP
    
    while N<maxiter and R_min2<r2<R_inf2:
        N+=1

        r2=sqrnorm(y[0:3])
        

        rkstep=fstep(r2)

        k1 = RK4f( y, L2)
        k2 = RK4f( y + 0.5*rkstep*k1, L2)
        k3 = RK4f( y + 0.5*rkstep*k2, L2)
        k4 = RK4f( y +     rkstep*k3, L2)

        increment = rkstep/6. * (k1 + 2*k2 + 2*k3 + k4)
                    
        y+= increment          
                
        X.append(y[0])
        Y.append(y[1])
        
    return X,Y


    

    
    


def render_1D():
    fig, ax = plt.subplots()
    ax.axis('equal')
    c1=plt.Circle((0, 0), 1., color='black',fill=False)
    c2=plt.Circle((0, 0), 1.5*R_schwarzschild, color='g',fill=False)

    adisk1=plt.Circle((0, 0), RAdisk_min, color='black',fill=False)
    adisk2=plt.Circle((0, 0), RAdisk_max, color='black',fill=False)


    ax.add_artist(c1)
    ax.add_artist(c2)
    ax.add_artist(adisk1)
    ax.add_artist(adisk2)

    T=np.linspace(0,1,10)
    for k in T:
        X,Y=sim_1D(np.array([-10,k,0]) ,np.array([np.pi/2,0.1]))
        plt.plot(X,Y,"*")

    plt.show()

if rendu=="1D":
    #start_time = time.time()
    render_1D()
    #print((time.time() - start_time))



def sim(point0,angle0,maxiter=10000):
    N=0
    #r=np.linalg.norm(pos)
    r2=sqrnorm(point0)
    oldz=0.

    y=np.zeros((6))
    
    y[0:3]=np.copy(point0)
    y[3:6]=asCartesian_u(angle0) #"vitesse"
   

    #vect=np.concatenate([point,velocity])
    L2=sqrnorm(np.cross(y[0:3],y[3:6])) #"moment cinétique"

    pixel_color=np.zeros((3))
    transparency_list=[]

    rkstep=STEP



    
    

    while N<maxiter and R_min2<r2<R_inf2:
        N+=1

        #rkstep=fstep(r2)


        oldz=y[2]

        r2=sqrnorm(y[0:3]) #Ne pas avoir à calculer de racine

        rkstep=fstep(r2)

        k1 = RK4f( y, L2)
        k2 = RK4f( y + 0.5*rkstep*k1, L2)
        k3 = RK4f( y + 0.5*rkstep*k2, L2)
        k4 = RK4f( y +     rkstep*k3, L2)

        increment = rkstep/6. * (k1 + 2*k2 + 2*k3 + k4)
                    
        y+= increment

        disk_crossing = np.logical_xor(oldz > 0., y[2] > 0.)#on traverse le plan z=0
        disk_distance = np.logical_and((r2 < RAdisk_max2), (r2 > RAdisk_min2)) #On l'a traversé la ou est le disque 
        disk = np.logical_and(disk_crossing,disk_distance)  

        if disk: #Il faut trouver l'intersection exacte avec le disque pour recuperer la texture
            lambdaa = - y[2]/y[5] #y[5] est la coordonné de la "vitesse" selon z
            coll_point=y[0:3]+lambdaa*y[3:6]
            r=np.linalg.norm(coll_point)
            if r<RAdisk_min:
                r=RAdisk_min
            if r>=RAdisk_max*0.999:
                r=RAdisk_max*0.999
            theta=np.arccos(coll_point[2]/r)
            a=int((r-RAdisk_min)*ADISK_X/(RAdisk_max-RAdisk_min))#pour ne pas depasser  
            b=int(ADISK_Y*np.mod((adisk_texture_rep*theta/np.pi+1)/2.,1)) #Pas besoin de s'occuper du modulo ici arctan2 dans -pi,pi
            transparency_list.append(adisk[a,b])# il faut cette liste pour pouvoir calculer la couleur du pixel si le rayon a croisé plusieurs fois le disque
                
                
        #vect=vect+step*eq(vect,L2) 

    n=len(transparency_list)
    
    angle=asSpherical_u(normalize(y[3:6]))
    if r2>R_min2:
        pixel_color=get_background_pixel(angle[0],angle[1]) #Le rayon  "a atteint l'infini" uniquement si la boucle precedente ne s'est pas arretee à cause de r<Rs
    if n>0:
        for i in range(n-1,-1,-1):
            alpha=np.max(transparency_list[i])#creer une transparence pour le disque (si la texture est noire il doit être totalement transparent), on peut quand meme voir le disque si le rayon a fini dans le trous noir
            pixel_color=alpha*transparency_list[i]+(1-alpha)*pixel_color 
    return pixel_color

    

def get_background_pixel(theta,phi):
    phi+=(1-np.sign(2*np.mod(theta-np.pi,2*np.pi)/(2*np.pi)-1))/2*np.pi
    a=np.abs(2*np.mod(theta-np.pi,2*np.pi)/(2*np.pi)-1)
    b=np.mod(np.pi-phi,2*np.pi)/(2*np.pi)
    return background[int(a*BACKGROUND_X),int(b*BACKGROUND_Y)] 


def render():
    start_time = time.time()
    res=np.zeros((HEIGHT,WIDTH,3))


    for i in range(HEIGHT):
        for j in range(WIDTH):

            #if HEIGHT*2.3/4<i<HEIGHT*2.6/4:      
            pos0=camera_pos
            directionvector=normalize([distance_viewplane,(WIDTH/2-j),(HEIGHT/2-i)])
            angle0=asSpherical_u(directionvector)+camera_dir

            res[i,j]=sim(pos0,angle0)
            
        print("%.2f pourcents \n " % (100*i/HEIGHT))
    print("Temps d execution : %s secondes " % (time.time() - start_time))
    return res


if rendu=="3D":

    res=render()

    plt.imsave("resultat.png",res)

    plt.imshow(res)
    plt.show()


#n=1
#T=np.linspace(-4*np.pi,4*np.pi,1000)
#Y=np.mod((n*T/np.pi+1)/2.,1)
#plt.plot(T,Y)
#plt.show()