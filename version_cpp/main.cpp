#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>
#include <string>
#include <fstream>
#include <vector>


using namespace std;

#include "raytracing.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define CHANNEL_NUM 3


#define MAXITER 6000

struct Rendu rdr;
struct Scene scn;




int background_width, background_height, background_bpp;
int adisk_width, adisk_height, adisk_bpp;

uint8_t* background = stbi_load("sphericalmap.png", &background_width, &background_height, &background_bpp, 3);
uint8_t* adisk = stbi_load("grille_adisk.png", &adisk_width, &adisk_height, &adisk_bpp, 3);




float sixth (float x)
{
    float tmp=x*x*x;
    tmp=tmp*tmp;
    return(tmp);
}

void asSpherical (float x, float y, float z,float* theta,float* phi)
{
    *theta=acos(z);
    *phi=atan2(y,x);
    
}

void asCartesian(float theta,float phi,float *x,float* y,float* z)
{
    float tmp=sin(theta);
    *x=tmp*cos(phi);
    *y=tmp*sin(phi);
    *z=cos(theta);
}


float mod(float x, float y)
{
   int resultat = floor (x / y);
   return x - (resultat * y);
}

float sqrnorm(float x,float y,float z){
    return x*x+y*y+z*z;
}
float norm(float x,float y,float z){
    return sqrt(x*x+y*y+z*z);
}

void normalise(float *x,float *y,float *z){
    float tmp=sqrt(sqrnorm(*x,*y,*z));
    *x/=tmp;
    *y/=tmp;
    *z/=tmp;
}


    


void sim(float x, float y,float z, float theta0,float phi0,int* pixelr,int* pixelg,int* pixelb)
{
    int N=0;
    float r2=x*x+y*y+z*z;
    float r6;
    float oldz; //#pour detecter le passage à travers le plan z=0 (pour tracer le disque)
    vector<float> pixel_transpr (rdr.maxtransparency);
    vector<float> pixel_transpg(rdr.maxtransparency);
    vector<float> pixel_transpb(rdr.maxtransparency);
    int k=0;

    float xp,yp,zp;

    asCartesian(theta0,phi0,&xp,&yp,&zp);   


    float rvectv_x = ( y * zp ) - ( yp * z );
	float rvectv_y = ( xp * z ) - ( x * zp );
	float rvectv_z = ( x * yp ) - ( xp * y );

    float L2=rvectv_x*rvectv_x+rvectv_y*rvectv_y+rvectv_z*rvectv_z;

    float tmp=-1.5*L2*rdr.step;
    bool  disk_crossing,disk_distance,disk;



    while ((N<MAXITER) && (rdr.R_min2<r2) && (r2<rdr.R_inf2))
    {
        N+=1;
        oldz=z;
        

        r2=sqrnorm(x,y,z);
        r6 = r2*r2*r2;

        
        x+=xp*rdr.step;
        y+=yp*rdr.step;
        z+=zp*rdr.step;

        

        xp+=  tmp * x / r6;
        yp+=  tmp * y / r6; 
        zp+=  tmp * z / r6;

        disk_crossing = (oldz>0) != (z > 0.) ;//on traverse le plan z=0
        disk_distance = (r2 < scn.RAdisk_max2) && (r2 > scn.RAdisk_min2); //On l'a traversé la ou est le disque 
        disk = disk_crossing && disk_distance;

        if (disk)
        { 
            float lambda = - z/zp; //#y[5] est la coordonné de la "vitesse" selon z
            float coll_x=x+lambda*xp;
            float coll_y=y+lambda*yp;
            float coll_z=z+lambda*zp;
            float r=norm(coll_x,coll_y,coll_z);

            disk=(r < scn.RAdisk_max) && (r > scn.RAdisk_min); //#reverification plus précise
            if (disk)
            {
                if (k<rdr.maxtransparency)
                {
                    k++;
                    float phi  =  atan2(coll_y/r,coll_x/r);
                    int cx=int((r-scn.RAdisk_min)*adisk_height/(scn.RAdisk_max-scn.RAdisk_min));
                    int cy=int(adisk_width*mod(scn.adisk_texture_rep*phi-M_PI,2.*M_PI)/(2.*M_PI));

                    int loc =(cx*adisk_width+cy)*CHANNEL_NUM;
                
                    pixel_transpr[k]==adisk[loc];
                    pixel_transpg[k]=adisk[loc+1];
                    pixel_transpb[k]=adisk[loc+2];
                }
                 
            }
        }

    }
    
    if (r2>rdr.R_min2){
        
        float theta;
        float phi;

        asSpherical (xp, yp, zp,&theta,&phi);


        if ((mod(theta-M_PI,2.*M_PI)/(M_PI)-1)<0.)
        {
        phi+=M_PI;
        }

        int cx=int(background_height*abs(2.*mod(theta-M_PI,2.*M_PI)/(2.*M_PI)-1.));
        int cy=int(background_width*mod(M_PI-phi,2.*M_PI)/(2.*M_PI));

        int loc=(cx*background_width+cy)*CHANNEL_NUM;

        *pixelr=background[loc];
        *pixelg=background[loc+1];
        *pixelb=background[loc+2];
    }

    for (int l = k; l >=0; --l)
    {
        float alpha=(float)pixel_transpg[l]/255.;

        *pixelr=int(alpha*(float)pixel_transpb[l]+(1-alpha)*(*pixelr));
        *pixelg=int(alpha*(float)pixel_transpb[l]+(1-alpha)*(*pixelg));
        *pixelb=int(alpha*(float)pixel_transpb[l]+(1-alpha)*(*pixelb)); 

        
    }        
        
}

void render(char* name){
    uint8_t* image = new uint8_t[rdr.width * rdr.height * CHANNEL_NUM];
    
    int index = 0;
    float x,y,z,xp0,yp0,zp0;
    float theta0,phi0;
    int pixelr,pixelg,pixelb;


     for (int j = 0; j <rdr.height; ++j)
     {
        for (int i = 0; i < rdr.width; ++i)
        {
        x=scn.camera.x;
        y=scn.camera.y;
        z=scn.camera.z;

        xp0=scn.distance_viewplane;
        yp0=((float)rdr.width/2.-i);
        zp0=((float)rdr.height/2.-j);

        normalise(&xp0,&yp0,&zp0);

        

        asSpherical(xp0,yp0,zp0,&theta0,&phi0);

        theta0+=scn.camera.theta;
        phi0+=scn.camera.phi;

        pixelr=0;
        pixelg=0;
        pixelb=0; 



        sim(x, y, z, theta0,phi0,&pixelr,&pixelg,&pixelb);

        //int coordx=int(a*height);
        //int coordy=int(b*width);

        //int loc=(j*background_width+i)*3;
        

        image[index++] = pixelr;
        image[index++] = pixelg;
        image[index++] = pixelb;

        }
        if (j%10==0){
        printf("file %s %.2f pourcents \n",name,100.*(float)j/((float)rdr.height));
        }
        
     }
  
    stbi_write_png(name, rdr.width, rdr.height, CHANNEL_NUM, image, rdr.width*CHANNEL_NUM);
}


    

    

int main() {

    std::string const monFichier("mesVariables.txt");
    ifstream monFlux(monFichier.c_str()) ;
    if (monFlux)
    {
    
    monFlux>>rdr.height;
    monFlux>>rdr.width;
    monFlux>>rdr.R_schwarzschild;
    monFlux>>rdr.R_inf; //distance à partir de laquelle on considere etre a l'infini
    monFlux>>rdr.step;
    monFlux>>rdr.maxtransparency;

    monFlux>>scn.RAdisk_min;
    monFlux>>scn.RAdisk_max;
    monFlux>>scn.FOV;
    
    monFlux>>scn.adisk_texture_rep;

    rdr.initRendering();//quleques calculs pour avoir les carrés de certaines qtité et le fov en radian etc..
    scn.initScene(rdr.width);

    //Rendre une image fixe:

    monFlux>>scn.camera.x;
    //scn.camera.y=0.;
    monFlux>>scn.camera.z;
    monFlux>>scn.camera.theta;
    monFlux>>scn.camera.phi;

    render("resultat.png");



    //rendre plusieurs images 

    /*float phicoordcamera=0;
    float thetacoordcamera=0;
    float rcoordcamera=13.;
    float xcamera,ycamera,zcamera;

    int nombreimage=150;

    char filename[10];


    for (int k=0;k<nombreimage;++k)
    {
        snprintf(filename, sizeof(filename), "%.3i%s", k+1, ".png");
        //printf ("%s \n",name);

        phicoordcamera=(float)k/((float)nombreimage)*2*M_PI*0.25+M_PI;
        thetacoordcamera=1.4937;
        asCartesian(thetacoordcamera,phicoordcamera,&xcamera,&ycamera,&zcamera);
        xcamera*=rcoordcamera;
        ycamera*=rcoordcamera;
        zcamera*=rcoordcamera;

        scn.camera.x=xcamera;
        scn.camera.y=ycamera;
        scn.camera.z=zcamera;
        scn.camera.theta=0.04;
        scn.camera.phi=(phicoordcamera-M_PI);

        render(filename);

    }*/

    stbi_image_free(background);
    stbi_image_free(adisk);

    return 0;
    }
    else
    {
        cout<<"ERREUR : Le fichier ne s'est pas ouvert correctement."<<endl;
    }
    
}






