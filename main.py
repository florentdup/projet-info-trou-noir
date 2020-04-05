import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import time

from PIL import Image
Image.MAX_IMAGE_PIXELS = 1000000000  

start_time = time.time()
background=mpimg.imread("eso0932a.tif")/255.
adisk=mpimg.imread("adisk.jpg")/255.
print("Temps de chargement image : %s secondes " % (time.time() - start_time))

BACKGROUND_X=background.shape[0]
BACKGROUND_Y=background.shape[1]

ADISK_X=adisk.shape[0]
ADISK_Y=adisk.shape[1]
alpha_adisk=0.95


HEIGHT=360*2 #x
WIDTH=640*2#y
FOV=50



RAdisk_min=2
RAdisk_max=5


R_schwarzschild=1 #inutile
R_min=1.*R_schwarzschild
R_inf=20
R_inf2=R_inf**2
R_inf6=R_inf**6
R_min6=1.


fov_rad=FOV*np.pi/180
distance_viewplane=WIDTH/(2*np.tan(fov_rad/2))


camera_pos=np.array([-10.,1.,0.5])
camera_dir=np.array([0.,0.05])#theta et phi

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
    return np.array([np.sin(thetaphi[0])*np.cos(thetaphi[1]),np.sin(thetaphi[0])*np.sin(thetaphi[1]),np.cos(thetaphi[0])])



def eq(x,L2):
        pos=x[:3]
        xp=x[3:]
        F=-3./2.*L2/sixth(pos)*pos
        return np.concatenate([xp,F])

def sim_1D(pos,angle,step=0.1,maxiter=10000):
    N=0
    r=np.linalg.norm(pos)
    velocity=asCartesian(angle)
    X=[]
    Y=[]

    vect=np.concatenate([pos,velocity])
    L2=sqrnorm(np.cross(pos,velocity))
    while N<maxiter and R_min/3<r<R_inf:
        r=np.linalg.norm(vect[:3])
        N+=1
        vect=vect+step*eq(vect,L2)
        X.append(vect[0])
        Y.append(vect[1])
        
    return X,Y


    

    
    


def render_1D():
    fig, ax = plt.subplots()
    ax.axis('equal')
    c1=plt.Circle((0, 0), 1., color='black',fill=False)
    c2=plt.Circle((0, 0), 1.5*R_schwarzschild, color='g',fill=False)
    ax.add_artist(c1)
    ax.add_artist(c2)

    T=np.linspace(0,1,50)
    for k in T:
        X,Y=sim_1D(np.array([-10,k,0]) ,np.array([1,np.pi/2,0.1]),0.01,10000)
        plt.plot(X,Y)

    plt.show()


#render_1D()


def sim(pos0,angle0,step=0.1,maxiter=10000):
    N=0
    #r=np.linalg.norm(pos)
    r6=sixth(pos0)
    pos=np.copy(pos0)
    velocity=asCartesian_u(angle0)

    #vect=np.concatenate([pos,velocity])
    L2=sqrnorm(np.cross(pos,velocity))

    pixel_color=np.zeros((3))
    transparency_list=[]
    

    while N<maxiter and R_min<r6<R_inf6:
        N+=1

        a=pos[2]

        r6=sixth(pos)
        pos+=velocity*step
        velocity+=-3./2.*L2/r6*pos*step

        

        if a*pos[2]<0:
            sph=asSpherical(pos)
            r=sph[0]
            if RAdisk_min<r<RAdisk_max:
                x=int((r-RAdisk_min)*ADISK_X/(RAdisk_max-RAdisk_min))
                y=int(ADISK_Y*(sph[2]/np.pi+1)/2.) #Pas besoin de s'occuper du modulo ici arctan2 dans -pi,pi
                transparency_list.append(adisk[x,y])# Pourvoir le disque si il est devant le trou noir
                
        #vect=vect+step*eq(vect,L2) 

    n=len(transparency_list)
    if r6>R_min6:
        angle=asSpherical_u(normalize(velocity))
        if n>0:
            transparency_list.append(get_background_pixel(angle[0],angle[1])) #A la fin de la liste 
            for i in range(n-1,-1,-1):
                pixel_color=alpha_adisk*transparency_list[i]+(1-alpha_adisk)*pixel_color 
        else:
            pixel_color=get_background_pixel(angle[0],angle[1])
    else:
        if  n>0:
            pixel_color=transparency_list[0] 
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
            pos=camera_pos
            directionvector=normalize([distance_viewplane,(WIDTH/2-j),(HEIGHT/2-i)])
            angle=asSpherical_u(directionvector)+camera_dir

            res[i,j]=sim(pos,angle,.1,10000)
            
        print("%.2f pourcents \n " % (100*i/HEIGHT))
    print("Temps d execution : %s secondes " % (time.time() - start_time))
    return res


res=render()

plt.imsave("resultat.png",res)

plt.imshow(res)
plt.show()



