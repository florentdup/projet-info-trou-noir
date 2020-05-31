#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>
#include <string.h>


using namespace std;

#include "fct.h"
#include "raytracing.h"


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define CHANNEL_NUM 3

#define MAXITER 6000



float FloatRand( float MaxVal );

void asCartesian(float theta,float phi,float *x,float* y,float* z);

void asSpherical (float x, float y, float z,float* theta,float* phi);

float mod(float x, float y);

float sqrnorm(float x,float y,float z);

float norm(float x,float y,float z);

void normalise(float* x,float* y,float* z);

float avg(float* array,int size);



struct Rendu rdr;
struct Scene scn;


float delta=.9;//"Taille angulaire (en pixel)" sur la direction des rayons d'un même groupe  




int adisk_width, adisk_height, adisk_bpp;
int colorScale_width, colorScale_height, colorScale_bpp;


uint8_t* adisk = stbi_load("adisk_upscaled.png", &adisk_width, &adisk_height, &adisk_bpp, 3);
uint8_t* colorScale = stbi_load("scale30000.png", &colorScale_width, &colorScale_height, &colorScale_bpp, 3); //Charger une echelle de couleur en fonction de la temperature

float Tmax=30000.;//Temperature max sur l'echelle en Kelvin
float Tmin=1000.;//Temperature min sur l'echelle en Kelvin


void normalisepixel(float* rgbR,float* rgbG,float* rgbB){
    float maxcomp=*rgbR;
    if (*rgbG>maxcomp){maxcomp=*rgbG;}
    if (*rgbB>maxcomp){maxcomp=*rgbB;}

    if (maxcomp>255.){
        *rgbR=255.*(*rgbR)/maxcomp;
        *rgbG=255.*(*rgbG)/maxcomp;
        *rgbB=255.*(*rgbB)/maxcomp;
    }  

}


void getRayDirection(float i,float j, float* xp ,float* yp,float* zp){

    *xp=scn.distance_viewplane;
    *yp=((float)rdr.width/2.-i);
    *zp=((float)rdr.height/2.-j);

    float theta0,phi0;

    normalise(xp,yp,zp);      
    asSpherical(*xp,*yp,*zp,&theta0,&phi0);  //Passer en corrodonnée spherique pour effectuer la rotation de la camera puis revenir en coord cartesiennes

    theta0+=scn.camera.theta;
    phi0+=scn.camera.phi;

    asCartesian(theta0,phi0,xp,yp,zp);
}


//Obtenir la couleur d'une étoile à partir de sa temperature (en corrigeant la luminosité)

//La variable brightness est un nombre quelconque pour amplifier la luminosité ou non
void getStarColor(float* rgbR,float* rgbG,float* rgbB,float temperature,float brightness){ 

    int c=int(((temperature-Tmin)*colorScale_width/(Tmax-Tmin)));

    int loc =(c)*CHANNEL_NUM;

    brightness*=1./(exp(29622.4/temperature)-1);//POur un corps noir dont le spectre est peu étalé

    *rgbR+=(float)colorScale[loc]*brightness;
    *rgbG+=(float)colorScale[loc+1]*brightness;
    *rgbB+=(float)colorScale[loc+2]*brightness;
      
}

//Même fonction mais cette fois sans la correction de luminosité
//La variable brightness est un nombre quelconque pour amplifier la luminosité ou non
void getDiskColor(float* rgbR,float* rgbG,float* rgbB,float temperature,float brightness){  

    int c=int(((temperature-Tmin)*colorScale_width/(Tmax-Tmin)));

    int loc =c*CHANNEL_NUM;

    *rgbR=(float)colorScale[loc]*brightness;
    *rgbG=(float)colorScale[loc+1]*brightness;
    *rgbB=(float)colorScale[loc+2]*brightness;
      
}

//Utilisé pour ajouter à l'image du fond celest le disque de poussière en superposant dans l'ordre les masques
void getFinalColor(float* rgbR,float* rgbG,float* rgbB,float* pixel_transpr,float* pixel_transpg,float* pixel_transpb,float* canalAlpha,int nbCollision){

     
    *rgbR*=200;//Il faut multiplier par un grand nombre à cause du 1./(exp(29622.4/temperature)-1) qui donne la luminosité
    *rgbG*=200;
    *rgbB*=200;
    normalisepixel(rgbR,rgbG,rgbB);

    for (int l = nbCollision-1; l >=0; --l)
    {
        pixel_transpr[l]*=70;
        pixel_transpg[l]*=70;
        pixel_transpb[l]*=70;//Il faut multiplier par un grand nombre à cause du 1./(exp(29622.4/temperature)-1) qui donne la luminosité

        normalisepixel(&pixel_transpr[l],&pixel_transpg[l],&pixel_transpb[l]);   

        *rgbR=canalAlpha[l]*pixel_transpr[l]+(1-canalAlpha[l])*(*rgbR);
        *rgbG=canalAlpha[l]*pixel_transpg[l]+(1-canalAlpha[l])*(*rgbG);
        *rgbB=canalAlpha[l]*pixel_transpb[l]+(1-canalAlpha[l])*(*rgbB);
        
    }
}
    

//Simulation complete d'un rayon (trajectoire+collision)
void sim(float x, float y,float z, float* xp,float* yp,float* zp,float* pixel_transpr,float* pixel_transpg,float* pixel_transpb,float* canalAlpha,bool* NotInBlackhole,int* nbCollision)
{
    int N=0;
    float r2=x*x+y*y+z*z;
    float r6;
    float oldz; //#pour detecter le passage à travers le plan z=0 (pour tracer le disque)


    int k=0;
 

    float rvectv_x = ( y * (*zp) ) - ( (*yp) * z );
	float rvectv_y = ( (*xp) * z ) - ( x * (*zp) );
	float rvectv_z = ( x * (*yp) ) - ( (*xp) * y );

    float L2=rvectv_x*rvectv_x+rvectv_y*rvectv_y+rvectv_z*rvectv_z;

    float tmp=-1.5*L2*rdr.step;
    bool  disk_crossing,disk_distance,disk;



    while ((N<MAXITER) && (rdr.R_min2<r2) && (r2<rdr.R_inf2))
    {
    
        N+=1;
        oldz=z;
        

        r2=sqrnorm(x,y,z);
        r6 = r2*r2*r2;

    
        x+=(*xp)*rdr.step;
        y+=(*yp)*rdr.step;
        z+=(*zp)*rdr.step;


        *xp+=  tmp * x / r6;
        *yp+=  tmp * y / r6; 
        *zp+=  tmp * z / r6;

        disk_crossing = (oldz>0) != (z > 0.) ;//on traverse le plan z=0
        disk_distance = (r2 < scn.RAdisk_max2) && (r2 > scn.RAdisk_min2); //On l'a traversé la ou est le disque 
        disk = disk_crossing && disk_distance;


        if (disk)
        { 
            float lambda = - z/(*zp); //#y[5] est la coordonné de la "vitesse" selon z
            float coll_x=x+lambda*(*xp);
            float coll_y=y+lambda*(*yp);
            float coll_z=z+lambda*(*zp);
            float r=norm(coll_x,coll_y,coll_z);

            disk=(r < scn.RAdisk_max) && (r > scn.RAdisk_min); //#reverification plus précise
            if (disk)
            {
                if (k<rdr.maxtransparency)
                {

                    float diskspeed_x,diskspeed_y,diskspeed_z;

                    
                    float phi  =  atan2(coll_y/r,coll_x/r);//Coordonné du point d'impact (en sphérique)

                    asCartesian(M_PI/2.,phi+M_PI/2.,&diskspeed_x,&diskspeed_y,&diskspeed_z); //Vitesse de rotation sens trigo,mvt circulaire
                    
                    float dirPhoton_x=*xp;
                    float dirPhoton_y=*yp;
                    float dirPhoton_z=*zp;

                    normalise(&dirPhoton_x,&dirPhoton_y,&dirPhoton_z);

                    float beta=scn.diskRotationSpeed(r);//Obtenir la vitesse des poussières en ce point
                    float costheta=(dirPhoton_x*diskspeed_x+dirPhoton_y*diskspeed_y+dirPhoton_z*diskspeed_z);//Obtenir le cosinus de l'angle entre le rayon et la vitesse des particules (Les vecteurs sont normés)

 
                    
                    //Temperature du disque à cet endroit
                    float Temp0=scn.diskTemp(r);
                    //Intensité du corps noir à cette temperature:
                    float luminosityfactor=1./(exp(29622.4/Temp0)-1);  //Le facteur d'intensité doit être calculé pour la temperature Temp0

                    //Décalage en fréquence environ égal à "décalage en temperature" (lambdamax*Temperature=cste)(uniquement pour la couleur, pas pour l'intensité)
                    float Temp=Temp0*(1-beta*costheta)/sqrt(1-beta*beta);//Effet doppler


                    Temp=Temp*sqrt((1.-1/r)); //Redshift gravitationnel simple à calculer ds la metrique de Schwarzschild (depend uniquement de la position d'emmision du rayon( sur le disque))
                    


                    float q=(Temp/Temp0);//Pour obtenir approximativement le rapport de fréquence nu/nu0 (frequence observée sur frequence d'emission)

                    luminosityfactor*=q*q*q*q; // (Intensité sur frequence au cube I/nu^3 est invariant, donc pour le spectre total décalage en puissance 4)

                    

                    getDiskColor(&pixel_transpr[k],&pixel_transpg[k],&pixel_transpb[k],Temp,luminosityfactor);//Obtenir la couleur


                    //Recuperer la texture du disque à cet endroit, qu'on utilise uniquemenet pour obtenir la transparence du disque
                    int cx=int((r-scn.RAdisk_min)*adisk_height/(scn.RAdisk_max-scn.RAdisk_min));
                    int cy=int(adisk_width*mod(scn.adisk_texture_rep*phi-M_PI,2.*M_PI)/(2.*M_PI));
                    int loc =(cx*adisk_width+cy)*CHANNEL_NUM;

                    canalAlpha[k]=(float)adisk[loc+1]/255.;;  //Choix: on utilise la composante verte (le +1) pour reconstituer la transparence
                                    
                    k++;
                }
                 
            }
        }

    }

    *NotInBlackhole=(r2>rdr.R_min2);
    *nbCollision=k;

    if(*NotInBlackhole){
        normalise(xp,yp,zp); //Ne pas effectuer inutilement cette opération
    }

}
    
//Simulation d'un rayon (trajectoire uniquement)
void sim_opt(float x, float y,float z, float* xp,float* yp,float* zp,bool* NotInBlackhole)
{
    int N=0;
    float r2=x*x+y*y+z*z;
    float r6;    
 


    float rvectv_x = ( y * (*zp) ) - ( (*yp) * z );
	float rvectv_y = ( (*xp) * z ) - ( x * (*zp) );
	float rvectv_z = ( x * (*yp) ) - ( (*xp) * y );

    float L2=rvectv_x*rvectv_x+rvectv_y*rvectv_y+rvectv_z*rvectv_z;

    float tmp=-1.5*L2*rdr.step;


    while ((N<MAXITER) && (rdr.R_min2<r2) && (r2<rdr.R_inf2))
    {
        N+=1;

        r2=sqrnorm(x,y,z);
        r6 = r2*r2*r2;

        
        x+=*xp*rdr.step;
        y+=*yp*rdr.step;
        z+=*zp*rdr.step;


        *xp+=  tmp * x / r6;
        *yp+=  tmp * y / r6; 
        *zp+=  tmp * z / r6;

    }

    *NotInBlackhole=(r2>rdr.R_min2);
    
    if(*NotInBlackhole){
        normalise(xp,yp,zp); //Ne pas effectuer inutilement cette opération
    }
}   

//Simulation d'un paquet de rayon
void sim_bundle(int i,int j,float x,float y,float z,float* pixelr,float* pixelg,float* pixelb){ //arguments: i et j position du pixel, x,y,z position du rayon dans l'espace
       
    float xpcentre0,ypcentre0,zpcentre0;
    float xpcentref,ypcentref,zpcentref;
    float xp0,yp0,zp0;
    float xpf,ypf,zpf;
    bool NotInBlackhole=true;
    int nbCollision;
    float maxdist2=0.;
    float maxdist2_0=0.;

    /*float mindist2=1000.;
    float mindist2_0=1000.;*/

    /*float Distances2[4];
    float Distances2_0[4];
    float dilatation=1.;*/
    


    float pixel_transpr[rdr.maxtransparency],pixel_transpg[rdr.maxtransparency],pixel_transpb[rdr.maxtransparency];
    float canalAlpha[rdr.maxtransparency];


    getRayDirection(i,j,&xpcentre0,&ypcentre0,&zpcentre0);  //Direction du rayon principal


    xpcentref=xpcentre0; 
    ypcentref=ypcentre0; 
    zpcentref=zpcentre0;

    sim(x,y,z,&xpcentref,&ypcentref,&zpcentref,pixel_transpr,pixel_transpg,pixel_transpb,canalAlpha,&NotInBlackhole,&nbCollision); //Simulation du rayon principal

    int Nray=0;

    //direction des rayons secondaires 
    //float theta=FloatRand(2.*M_PI);    
    
    //M_PI/80.*(float)((i*100/rdr.width)+(j*100/rdr.height));
    
    //float c=cos(*angle);
    //float s=sin(*angle);


    while (Nray<4 && NotInBlackhole){

        if (Nray==0){ getRayDirection(i+delta,j,&xp0,&yp0,&zp0); }
        if (Nray==1){ getRayDirection(i-delta,j,&xp0,&yp0,&zp0); }
        if (Nray==2){ getRayDirection(i,j+delta,&xp0,&yp0,&zp0); }
        if (Nray==3){ getRayDirection(i,j-delta,&xp0,&yp0,&zp0); }

        //float noiseX=0.;
        //float noiseY=0.;

        
        //Version bruit sur les directions
        /*while (noiseX*noiseX+noiseY*noiseY<0.4){
            noiseX=FloatRand(2.*noisedelta)-noisedelta;
            noiseY=FloatRand(2.*noisedelta)-noisedelta;
        }*/

        //Version ou la croix tourne d'un angle theta potentiellement pas le même pour tous les pixels
        //getRayDirection(i+delta+noiseX,j+noiseY,&xp0,&yp0,&zp0);
        /*if (Nray==0){ getRayDirection(i+delta*c,j+delta*s,&xp0,&yp0,&zp0); }
        if (Nray==1){ getRayDirection(i-delta*s,j+delta*c,&xp0,&yp0,&zp0); }
        if (Nray==2){ getRayDirection(i-delta*c,j-delta*s,&xp0,&yp0,&zp0); }
        if (Nray==3){ getRayDirection(i+delta*s,j-delta*c,&xp0,&yp0,&zp0); }*/
        Nray++;


        xpf=xp0;
        ypf=yp0;
        zpf=zp0;

        sim_opt(x,y,z,&xpf,&ypf,&zpf,&NotInBlackhole); //Simulation des rayons secondaires

        if(NotInBlackhole){
            //Calcul de la taille de la zone d'impact des differents rayons (zone en terme de direction)
            //On trouve le rayon maxdist du cercle centré sur la direction du rayon principal, qui contient les autres rayons (tjrs en terme de direction)
            //Calcul en coord cartesienne pour éviter les discontinuités liées au modulo 2PI

            float dx0=xpcentre0-xp0; 
            float dy0=ypcentre0-yp0; 
            float dz0=zpcentre0-zp0;

            float dist2_0=sqrnorm(dx0,dy0,dz0);


            float dx=xpcentref-xpf;  
            float dy=ypcentref-ypf;
            float dz=zpcentref-zpf;

            float dist2=sqrnorm(dx,dy,dz);

            if (maxdist2_0<dist2_0){ maxdist2_0=dist2_0; }
            if (maxdist2<dist2){ maxdist2=dist2; }

            /*if (Nray==0){mindist2=dist2;}
            else if (mindist2>dist2){ mindist2=dist2; }*/

            /*if (Nray==0){mindist2_0=dist2_0;}
            else if (mindist2_0>dist2_0){ mindist2_0=dist2_0; }

            if (mindist2_0>dist2_0){ mindist2_0=dist2_0; }*/

            

            /*Distances2_0[Nray]=dist2_0;
            Distances2[Nray]=dist2;*/

            
            }
        }


        //J'ai essayé d'autres méthodes pour estimer la dispertion du faiscaux

        /*for (int l=0;l<4;++l){
                if (Distances2[l]>0){dilatation+=sqrt(Distances2_0[l]/Distances2[l]);}
        }
        dilatation/=4.;

        float R2Init=avg(Distances2_0,4);
        float R2f=R2Init*dilatation*dilatation;*/



        if(NotInBlackhole){//Si pas un seul des rayons du groupe n'est tombé dans le trou noir

            float maxdist=sqrt(maxdist2);  
            float sqrdistanceToStar;//sqrnorm(scn.starx[i]-xpcentref,scn.stary[i]-ypcentref,scn.starz[i]-zpcentref);

            float c1,c2,c3;

            
            for(int i=0;i<Nstar;++i){//Tres couteux 

                c1=abs(scn.starx[i]-xpcentref);

                if (c1<maxdist){
                    c2=abs(scn.stary[i]-ypcentref);
                    if (c2<maxdist){
                        c3=abs(scn.starz[i]-zpcentref);
                        if (c3<maxdist){
                            sqrdistanceToStar=sqrnorm(c1,c2,c3);
                            if (sqrdistanceToStar<maxdist2){ 
                                getStarColor(pixelr,pixelg,pixelb,scn.starTemp[i],scn.starBrightness[i]); //(Récupere la couleur et effectue l'addition composante par composante)
    
                                
                            }

                        }
                    }
                }                               
                  
            }

            /*float I[4];

            for (int l=0;l<4;++l){
                if (Distances2[l]>0){I[l]=sqrt(Distances2_0[l]/Distances2[l]);}
            }

            float minL=I[0];
            float maxL=I[0];
            float avgL=I[0];
            for (int l=1;l<4;++l){
                if (minL>I[l]){minL=I[l];}
                if (maxL<I[l]){maxL=I[l];}
                avgL+=I[l];
            }
            avgL/=4.;*/

            
            //float b=dilatation*dilatation;
            float b=maxdist2_0/maxdist2;
            *pixelr*=b;
            *pixelg*=b;
            *pixelb*=b;
        }

        getFinalColor(pixelr,pixelg,pixelb,pixel_transpr,pixel_transpg,pixel_transpb,canalAlpha,nbCollision);// Ajouter les superpositions liées au disque de poussière (calcul de transparence)
       
}
        


void render(const char* name){
    float* image = new float[rdr.width * rdr.height * CHANNEL_NUM]; //Tableau qui va contenir les données
    
    int index = 0;
    float x,y,z;
    float pixelr,pixelg,pixelb;

    //int maxintensity=0;

    float angle=0.;


     for (int j = 0; j <rdr.height; ++j)
     {
        for (int i = 0; i < rdr.width; ++i)
        {    

        //Position initial du rayon à l'itération 0
        x=scn.camera.x;
        y=scn.camera.y;
        z=scn.camera.z; 

        pixelr=0.;
        pixelg=0.;
        pixelb=0.; 

        sim_bundle(i,j,x,y,z,&pixelr,&pixelg,&pixelb);


        image[index++] = pixelr;
        image[index++] = pixelg;
        image[index++] = pixelb;

        }
        if (j%10==0){
        printf("file %s %.2f pourcents \n",name,100.*(float)j/((float)rdr.height));  
        }   
     }


     //Postprocess possible ici

     index=0;
     
     uint8_t* image_byte = new uint8_t[rdr.width * rdr.height * CHANNEL_NUM]; //Image finale

     

     for (int j = 0; j <rdr.height; ++j)
     {
        for (int i = 0; i < rdr.width; ++i)
        {

            image_byte[index]=char(image[index]);//Puis on place la coordonne dans l'image finale
            index++;
            
            image_byte[index]=char(image[index]);
            index++;
         
            image_byte[index]=char(image[index]);
            index++;
            
        }
    }
  
    stbi_write_png(name, rdr.width, rdr.height, CHANNEL_NUM, image_byte, rdr.width*CHANNEL_NUM);
}

    

int main() {

    rdr.height=1080;
    rdr.width=1920;
    rdr.R_schwarzschild=1.;
    rdr.R_inf=15.; //distance à partir de laquelle on considere etre a l'infini
    rdr.step=0.05;//Pas d'intégration
    rdr.maxtransparency=4;//Nombre de collision max pour un rayon

    scn.RAdisk_min=1.5;
    scn.RAdisk_max=10.;
    scn.betamax=0.3;//Vitesse du disque au plus proche du trou noir (c'est la qu'est la vitesse max pour un profil de vitesse en /r^-1/2)
    scn.TdiskMax=15000.;//temperature du disque au plus proche du trou noir

    scn.FOV=40.;
    
    scn.adisk_texture_rep=4.; //Répéter la texture du dique (en longeur pour ne pas qu'elle soit pixelisée)

    rdr.initRendering();//quelques calculs pour avoir les carrés de certaines qtité et le fov en radian etc..
    scn.initScene(rdr.width);

    scn.generateSky();

    //Rendre une image fixe:
    scn.camera.x=-13.;
    scn.camera.y=0.;
    scn.camera.z=1.;
    scn.camera.theta=0.04;
    scn.camera.phi=0.;

   


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

    stbi_image_free(adisk);
    stbi_image_free(colorScale);

    return 0;
}