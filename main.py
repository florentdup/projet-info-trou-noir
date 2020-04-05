import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import time

from PIL import Image
Image.MAX_IMAGE_PIXELS = 1000000000  

start_time = time.time()
background=mpimg.imread("eso_space.jpg")/255.
print("Temps de chargement image : %s secondes " % (time.time() - start_time))

BACKGROUND_X=background.shape[0]
BACKGROUND_Y=background.shape[1]


HEIGHT=360*2 #x
WIDTH=640*2 #y
FOV=80

SIZE=1.


R_schwarzschild=1#inutile
Rinf=20


fov_rad=FOV*np.pi/180
distance_viewplane=WIDTH*SIZE/(2*np.tan(fov_rad/2))


camera_pos=np.array([-10.,0.,0.])
camera_dir=np.array([0.,0.])#theta et phi

def normalize(x):
    return x/np.linalg.norm(x)

def sqrnorm(x):
    return np.dot(x,x)

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


def eq(x,L2):
        pos=x[:3]
        xp=x[3:]
        r=np.linalg.norm(pos)
        F=-3./2.*L2/r**6*pos
        return np.concatenate([xp,F])

def sim_1D(pos,angle,step=0.1,maxiter=10000):
    N=0
    r=np.linalg.norm(pos)
    velocity=asCartesian(angle)
    X=[]
    Y=[]

    vect=np.concatenate([pos,velocity])
    L2=sqrnorm(np.cross(pos,velocity))
    while N<maxiter and 1.5*R_schwarzschild<r<Rinf:
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
        X,Y=sim_1D(np.array([-10,k,0]) ,np.array([1,np.pi/2,0.1]),0.1,1000)
        plt.plot(X,Y)

    plt.show()


#render_1D()


R_max=1.5*R_schwarzschild

def sim(pos,angle,step=0.1,maxiter=10000):
    N=0
    r=np.linalg.norm(pos)
    velocity=asCartesian(angle)

    vect=np.concatenate([pos,velocity])
    L2=sqrnorm(np.cross(pos,velocity))

    while N<maxiter and R_max<r<Rinf:
        r=np.linalg.norm(vect[:3])
        N+=1
        vect=vect+step*eq(vect,L2) 
    return asSpherical(vect[3:]),r

    

def getpixel(theta,phi):
    phi+=(1-np.sign(2*np.mod(theta-np.pi,2*np.pi)/(2*np.pi)-1))/2*np.pi
    a=np.abs(2*np.mod(theta-np.pi,2*np.pi)/(2*np.pi)-1)
    b=np.mod(np.pi-phi,2*np.pi)/(2*np.pi)
    return background[int(a*BACKGROUND_X),int(b*BACKGROUND_Y)] 


def render():
    start_time = time.time()
    res=np.zeros((HEIGHT,WIDTH,3))


    for i in range(HEIGHT):
        for j in range(WIDTH):       
            pos=np.copy(camera_pos)
            directionvector=normalize([distance_viewplane,(WIDTH/2-j)*SIZE,(HEIGHT/2-i)*SIZE])
            angle=(asSpherical(directionvector)[1:])+camera_dir

            angle_2=np.array([1,angle[0],angle[1]])

            angle_final,ray_distance=sim(pos,angle_2,0.1,1000)
            theta=angle_final[1]
            phi=angle_final[2]
            if ray_distance>R_max:
                res[i,j]=getpixel(theta,phi)#sinon le pixel reste noir
            
        print("%.2f pourcents \n " % (100*i/HEIGHT))
    print("Temps d execution : %s secondes " % (time.time() - start_time))
    return res


res=render()

plt.imsave("resultat.png",res)

plt.imshow(res)
plt.show()





