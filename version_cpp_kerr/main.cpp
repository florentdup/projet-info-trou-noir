#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>
#include <string.h>


using namespace std;

#include "fct.h"
#include "raytracing.h"
#include "par_for.h"


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define CHANNEL_NUM 3

#define MAXITER 500000
static float HM_PI=M_PI/2.;




float FloatRand( float MaxVal );

void sphericalToCartesian(float theta,float phi,float *x,float* y,float* z);

void cartesianToSpherical (float x, float y, float z,float* theta,float* phi);

float mod(float x, float y);

float sqrnorm(float x,float y,float z);

float norm(float x,float y,float z);

void normalise(float* x,float* y,float* z);

float avg(float* array,int size);

static const int n=6;





struct Rendu rdr;
struct Scene scn;
struct Blackhole bh;
struct Disk disk;







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
    cartesianToSpherical(*xp,*yp,*zp,&theta0,&phi0);  //Passer en corrodonnée spherique pour effectuer la rotation de la camera puis revenir en coord cartesiennes

    theta0+=scn.camera.theta;
    phi0+=scn.camera.phi;

    sphericalToCartesian(theta0,phi0,xp,yp,zp);
}


//Obtenir la couleur d'une étoile à partir de sa temperature (en corrigeant la luminosité)

//La variable brightness est un nombre quelconque pour amplifier la luminosité ou non
void getStarColor(float* rgbR,float* rgbG,float* rgbB,float temperature,float brightness){ 

    int c=int(((temperature-Tmin)*colorScale_width/(Tmax-Tmin)));

    int loc =c*CHANNEL_NUM;

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

    brightness*=1./(exp(29622.4/temperature)-1);//POur un corps noir dont le spectre est peu étalé

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


void cartesianToBl(float x,float y,float z,float* r,float* theta,float* phi){//Faux

    float r2=x*x+y*y+z*z;
    *r=sqrt(r2);

    *phi=atan2(y,x);
    *theta=acos(z/(*r));
 
}

void BlToCartesian(float r,float theta,float phi,float* x,float*y,float*z){

    float sintheta=sin(theta);
    float costheta=cos(theta);
    float cosphi=cos(phi);
    float sinphi=sin(phi);
    float temp=sintheta*sqrt(r*r+bh.a2);

    *x=temp*cosphi;
    *y=temp*sinphi;
    *z=r*costheta;

}

void getBlDirection(float x,float y,float z,float xp,float yp,float zp,float* rdot,float* thetadot,float* phidot){//xp,yp,zp normalisé

    float r0,theta0,phi0;
    float r1,theta1,phi1;
    
    cartesianToBl(x,y,z,&r0,&theta0,&phi0);
    cartesianToBl(x+xp*rdr.defaultstep,y+yp*rdr.defaultstep,z+zp*rdr.defaultstep,&r1,&theta1,&phi1);
    float dr=r1-r0;
    float dtheta=theta1-theta0;
    float dphi=phi1-phi0;

    *rdot=dr/rdr.defaultstep;
    *thetadot=dtheta/rdr.defaultstep;
    *phidot=dphi/rdr.defaultstep;

}

/* Coupled differential equations describing motion of photon */
void geodesic(float L,float kappa,float *y, float *dydx)
{
	float r, theta, pr, ptheta;

	r = y[0];
	theta = y[1];
	pr = y[4];
	ptheta = y[5];

	float r2 = r*r;
	float twor = 2.0*r;

	float sintheta, costheta;
	sintheta=sin(theta);
    costheta=cos(theta);

    //*side=(costheta>0.);
	
	float cos2 = costheta*costheta;
	float sin2 = sintheta*sintheta;

	float sigma = r2+bh.a2*cos2;
	float delta = r2-twor+bh.a2;
	float sd = sigma*delta;
	float siginv = 1.0/sigma;
	float bot = 1.0/sd;

	/* Prevent problems with the axis */
	if (sintheta < 1e-8)
	{
		sintheta = 1e-8;
		sin2 = 1e-16;
	}
    

	dydx[0] = -pr*delta*siginv;
	dydx[1] = -ptheta*siginv;
	dydx[2] = -(twor*bh.a+(sigma-twor)*L/sin2)*bot;
	dydx[3] = -(1.0+(twor*(r2+bh.a2)-twor*bh.a*L)*bot);
	dydx[4] = -(((r-1.0)*(-kappa)+twor*(r2+bh.a2)-2.0*bh.a*L)*bot-2.0*pr*pr*(r-1.0)*siginv);
	dydx[5] = -sintheta*costheta*(L*L/(sin2*sin2)-bh.a2)*siginv;
}

/* Initial Conditions for Ray */
void initial(float r0,float theta0,float* L,float* kappa,float *y0, float *ydot0, float x, float y)
{
	y0[0] = r0;
	y0[1] = theta0;
	y0[2] = 0;
	y0[3] = 0;
	y0[4] = cos(y)*cos(x);
	y0[5] = sin(y)/r0;

	float sintheta, costheta;
	sintheta=sin(theta0);
    costheta=cos(theta0);
	float cos2 = costheta*costheta;
	float sin2 = sintheta*sintheta;

	float rdot0 = y0[4];
	float thetadot0 = y0[5];

	float r2 = r0 * r0;
	float sigma = r2 + bh.a2*cos2;
	float delta = r2 - 2.0 * r0 + bh.a2;
	float s1 = sigma - 2.0 * r0;

	y0[4]= rdot0*sigma/delta;
	y0[5]= thetadot0*sigma;

	ydot0[0] = rdot0;
	ydot0[1] = thetadot0;
	ydot0[2] = cos(y)*sin(x)/(r0*sin(theta0));

	float phidot0 = ydot0[2];
	float energy2 = s1*(rdot0*rdot0/delta+thetadot0*thetadot0)
					+ delta*sin2*phidot0*phidot0;

	float energy = sqrt(energy2);

	/* Rescale */
	y0[4] = y0[4]/energy;
	y0[5] = y0[5]/energy;

	/* Angular Momentum with E = 1 */
	*L = ((sigma*delta*phidot0-2.0*bh.a*r0*energy)*sin2/s1)/energy;

	*kappa = y0[5]*y0[5]+bh.a2*sin2+(*L)*(*L)/sin2;

	/* Hack - make sure everything is normalized correctly by a call to geodesic */

 
	geodesic(*L,*kappa,y0, ydot0);
}




void getPosAndDir(float* y,float* dydx,float* xpos,float* ypos,float* zpos,float* xp,float* yp,float* zp){
    
    BlToCartesian(y[0],y[1],y[2],xpos,ypos,zpos);

    float costheta=cos(y[1]);
    float sintheta=sin(y[1]);
    float cosphi=cos(y[2]);
    float sinphi=sin(y[2]);
    float r2=y[0]*y[0];
    float R2=r2+bh.a2;
    float R=sqrt(r2+bh.a2);
            
    *xp=R*dydx[1]*cosphi*costheta-R*dydx[2]*sinphi*sintheta+sintheta*cosphi*y[0]*dydx[0]/R;
    *yp=R*dydx[1]*sinphi*costheta+R*dydx[2]*cosphi*sintheta+sintheta*sinphi*y[0]*dydx[0]/R;
    *zp=-y[0]*sintheta*dydx[1]+dydx[0]*costheta;
    normalise(xp,yp,zp);
}

void getDir(float* y,float* dydx,float* xp,float* yp,float* zp){
    
    float costheta=cos(y[1]);
    float sintheta=sin(y[1]);
    float cosphi=cos(y[2]);
    float sinphi=sin(y[2]);
    float r2=y[0]*y[0];
    float R2=r2+bh.a2;
    float R=sqrt(r2+bh.a2);
            
    *xp=R*dydx[1]*cosphi*costheta-R*dydx[2]*sinphi*sintheta+sintheta*cosphi*y[0]*dydx[0]/R;
    *yp=R*dydx[1]*sinphi*costheta+R*dydx[2]*cosphi*sintheta+sintheta*sinphi*y[0]*dydx[0]/R;
    *zp=-y[0]*sintheta*dydx[1]+dydx[0]*costheta;
    normalise(xp,yp,zp);
}


    

//Simulation complete d'un rayon (trajectoire+collision)
void sim(float posx, float posy,float posz, float xpixel,float ypixel,float* xp,float* yp,float* zp,float* pixel_transpr,float* pixel_transpg,float* pixel_transpb,float* canalAlpha,bool* NotInBlackhole,int* nbCollision)
{
    int N=0;
    int k=0;    

    
    
    

    float oldz; //#pour detecter le passage à travers le plan z=0 (pour tracer le disque)

    bool  zSignChange,diskDistance,diskCollision;

    //CI et constantes:
    float r0,theta0,phi0;

    cartesianToBl(posx,posy,posz,&r0,&theta0,&phi0);



    
    float y[6];
    float oldy[3];
    float dydx[6];

    /*float ak[6];
    float y1[6];
    float y2[6];
    float y3[6];
    float y4[6];
    float ytemp[6];*/
    

    float L,kappa;

    initial(r0,theta0,&L,&kappa,y,dydx,xpixel,ypixel);

    //bool oldside=true;
    //bool side=true;




    /*float rdot0,phidot0,thetadot0;

    //getBlDirection(x,y,z,*xp,*yp,*zp,&rdot0,&thetadot0,&phidot0);//obtenir la derivée des coord du rayon 

    rdot0=cos(ypixel)*cos(xpixel);
    phidot0=cos(ypixel)*sin(xpixel)/(r0*sin(theta0));
    thetadot0=sin(ypixel)/r0;

    //Constantes:

    float sintheta0=sin(theta0);
    float costheta0=cos(theta0);

    float sintheta0_2=sintheta0*sintheta0;
    float costheta0_2=costheta0*costheta0;


    float delta0=r0*r0-2.*r0+bh.a2;
    float sigma0=r0*r0+bh.a2*costheta0_2;


    float C=sigma0-2.*r0;
    float energy2 = C*(rdot0*rdot0/delta0+thetadot0*thetadot0)+ delta0*sintheta0_2*phidot0*phidot0;

	float energy = sqrt(energy2);

    float thetadot0sigma=thetadot0*sigma0;

    float L = ((sigma0*delta*phidot0-2.0*bh.a*r0*energy)*sintheta0_2/C)/energy;

    float L2=L*L;

    float K=thetadot0sigma*thetadot0sigma+bh.a2*sintheta0_2+L2/sintheta0_2;



    float delta=delta0;
    float sigma=sigma0;
    float costheta=costheta0;
    float sintheta=sintheta0;
    float costheta2=costheta0_2;
    float sintheta2=sintheta0_2;
    float sigmaInv=1./sigma0;
    float deltasigma=delta*sigma;

    float r=r0;
    float r2=r*r;
    float twor=2.*r;

    float phi=phi0;
    float theta=theta0;
    float pr=rdot0*sigma/(delta*energy);
    float ptheta=sigma/energy;


    //printf("%f %f %f %f \n",r,energy,sigma,delta);*/


    
    //float r2=y[0]*y[0];


    float xpos,ypos,zpos;


    float currentStep=rdr.step(y[0]);

    int l=0;;   

    


    while ((N<MAXITER) && (rdr.R_min<y[0]) && (y[0]<rdr.R_inf))
    {
    
        N+=1;

        //oldside=side;

        oldy[0]=y[0];
        oldy[1]=y[1];
        oldy[2]=y[2];

        
        
        currentStep=rdr.step(y[0]);



        for (int l = 0; l < 6; l++)
	    {
		    float hdydx = currentStep * dydx[l];
		    y[l] = y[l] + hdydx;
	    }

	    geodesic(L,kappa,y, dydx);

        /*geodesic(L,kappa,y, ak);
        for (l=0;l<6;++l)
        {
            y1[l]=ak[l];
            ytemp[l]=y[l]+.5*currentStep*ak[l];
        }
        geodesic(L,kappa,ytemp,ak);
        for (l=0;l<6;++l)
        {
            y2[l]=2.*ak[l];
            ytemp[l]=y[l]+.5*currentStep*ak[l];
        }
        geodesic(L,kappa,ytemp,ak);
        for (l=0;l<6;++l)
        {
            y3[l]=2.*ak[l];
            ytemp[l]=y[l]+currentStep*ak[l];
        }
        geodesic(L,kappa,ytemp,ak);
        for (l=0;l<6;++l)
        {
            y4[l]=ak[l];
            y[l]=y[l]+currentStep/6*(y1[l]+2.*y2[l]+2.*y3[l]+y4[l]);
        }*/
    

        /*r2=r*r;
        twor=2.0*r;

        costheta=cos(theta);
        sintheta=sin(theta);
        costheta2=costheta*costheta;
        sintheta2=sintheta*sintheta;
        
        delta=r2-twor+bh.a2;
        sigma=r2+bh.a2*costheta2;
        sigmaInv=1./sigma;

        deltasigma=delta*sigma;

        r=r+delta*sigmaInv*pr*ds;
        theta=theta+sigmaInv*ptheta*ds;

        phi=phi+(twor*bh.a+(sigma-twor)*L/sintheta2)/deltasigma*ds;

        pr=pr+((twor*(r2+bh.a2)-(r-1.)*K-2.*bh.a*L)/deltasigma-pr*pr*(r-1.)/sigma)*ds;

        ptheta=ptheta+sintheta*costheta*L2/bh.a2-bh.a2/sigma*ds;*/


        

        /*r2=sqrnorm(x,y,z);
        r6 = r2*r2*r2;

    
        x+=(*xp)*rdr.step;
        y+=(*yp)*rdr.step;
        z+=(*zp)*rdr.step;


        *xp+=  tmp * x / r6;
        *yp+=  tmp * y / r6; 
        *zp+=  tmp * z / r6;*/







        zSignChange = (oldy[1]>HM_PI) != (y[1]>HM_PI) ;//on traverse le plan z=0
        diskDistance = (y[0] < disk.R_max) && (y[0] > disk.R_min); //On l'a traversé la ou est le disque 
        diskCollision = zSignChange && diskDistance;


        if (diskCollision)
        { 
            getPosAndDir(y,dydx,&xpos,&ypos,&zpos,xp,yp,zp); 
            //BlToCartesian(y[0],y[1],y[2],&xpos,&ypos,&zpos);
            //zp=1.;
            /*float h=norm(xp,yp,zp);
            if (h<1e-5){
                zp+=0.001;
            }*/
            
            float lambda = - zpos/(*zp); 
            float coll_x=xpos+lambda*(*xp);
            float coll_y=ypos+lambda*(*yp);
            float coll_z=zpos+lambda*(*zp);
            float r=norm(coll_x,coll_y,coll_z);

            diskCollision=(r < disk.R_max) && (r > disk.R_min); //#reverification plus précise
            if (diskCollision)
            {
                if (k<rdr.maxtransparency)
                {

                    float diskspeed_x,diskspeed_y,diskspeed_z;

                    
                    float phi  =  atan2(coll_y/r,coll_x/r);//Coordonné du point d'impact (en sphérique)

                    sphericalToCartesian(M_PI/2.,phi+M_PI/2.,&diskspeed_x,&diskspeed_y,&diskspeed_z); //Vitesse de rotation sens trigo,mvt circulaire
                    
                    float dirPhoton_x=*xp;
                    float dirPhoton_y=*yp;
                    float dirPhoton_z=*zp;

                    normalise(&dirPhoton_x,&dirPhoton_y,&dirPhoton_z);

                    float beta=disk.RotationSpeed(r);//Obtenir la vitesse des poussières en ce point
                    float costhetadoppler=-(dirPhoton_x*diskspeed_x+dirPhoton_y*diskspeed_y+dirPhoton_z*diskspeed_z);//Obtenir le cosinus de l'angle entre le rayon et la vitesse des particules (Les vecteurs sont normés)

 
                    
                    //Temperature du disque à cet endroit
                    float Temp0=disk.Temp(r);
                    //Intensité du corps noir à cette temperature:
                    //float luminosityfactor=1./(exp(29622.4/Temp0)-1);  //Le facteur d'intensité doit être calculé pour la temperature Temp0

                    //Décalage en fréquence environ égal à "décalage en temperature" (lambdamax*Temperature=cste)(uniquement pour la couleur, pas pour l'intensité)
                    float Temp=Temp0*(1+beta*costhetadoppler)/sqrt(1-beta*beta);//Effet doppler
                    //float Temp=Temp0;


                    //Temp=Temp*sqrt((1.-1/r)); //Redshift gravitationnel simple à calculer ds la metrique de Schwarzschild (depend uniquement de la position d'emmision du rayon( sur le disque))
                    
                    

                    //float q=(Temp/Temp0);//Pour obtenir approximativement le rapport de fréquence nu/nu0 (frequence observée sur frequence d'emission)
                    //luminosityfactor*=q*q*q*q; // (Intensité sur frequence au cube I/nu^3 est invariant, donc pour le spectre total décalage en puissance 4)


                    //Recuperer la texture du disque à cet endroit, qu'on utilise uniquemenet pour obtenir la transparence du disque
                    int cx=int((r-disk.R_min)*adisk_height/(disk.R_max-disk.R_min));
                    int cy=int(adisk_width*mod(disk.texture_rep*phi-M_PI,2.*M_PI)/(2.*M_PI));
                    int loc =(cx*adisk_width+cy)*CHANNEL_NUM;

                    getDiskColor(&pixel_transpr[k],&pixel_transpg[k],&pixel_transpb[k],Temp,1.);//Obtenir la couleur

                    canalAlpha[k]=(float)adisk[loc+1]/255.;  //Choix: on utilise la composante verte (le +1) pour reconstituer la transparence
                                    
                    k++;
                }
                 
            }
        }

    }
    //printf("%d \n",N);

    *NotInBlackhole=(y[0]>rdr.R_min);
    *nbCollision=k;

    if(*NotInBlackhole){
        getDir(y,dydx,xp,yp,zp);
        normalise(xp,yp,zp); //Ne pas effectuer inutilement cette opération
    }

}
    
//Simulation d'un rayon (trajectoire uniquement)
void sim_opt(float posx, float posy,float posz, float xpixel,float ypixel,float* xp,float* yp,float* zp,bool* NotInBlackhole)
{
    int N=0;
    
    float r0,theta0,phi0;

    cartesianToBl(posx,posy,posz,&r0,&theta0,&phi0);


    
    float y[6];
    float dydx[6];

    float L,kappa;

    initial(r0,theta0,&L,&kappa,y,dydx,xpixel,ypixel);


    float currentStep=rdr.step(y[0]);

    int l=0;  


    while ((N<MAXITER) && (rdr.R_min<y[0]) && (y[0]<rdr.R_inf))
    {
        N+=1;

        currentStep=rdr.step(y[0]);

        for (int l = 0; l < 6; l++)
	    {
		    float hdydx = currentStep * dydx[l];
		    y[l] = y[l] + hdydx;
	    }

	    geodesic(L,kappa,y, dydx);


    }

    *NotInBlackhole=(y[0]>rdr.R_min);


    
    
    if(*NotInBlackhole){
        getDir(y,dydx,xp,yp,zp);
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


    //getRayDirection(i,j,&xpcentre0,&ypcentre0,&zpcentre0);  //Direction du rayon principal


    xpcentref=xpcentre0; 
    ypcentref=ypcentre0; 
    zpcentref=zpcentre0;

  

    float range = 0.08 * 20. / (rdr.width - 1.0);




    float xpixel=-(i - (rdr.width + 1.0) / 2)*range;  
    float ypixel=-(j - (rdr.height + 1.0) / 2)*range;

    //Eviter de passer par les poles
    if ((xpixel<0) && (xpixel>.001)){xpixel=-.001;}
    if ((xpixel>0) && (xpixel<.001)){xpixel=.001;}



    sim(x,y,z,xpixel,ypixel,&xpcentref,&ypcentref,&zpcentref,pixel_transpr,pixel_transpg,pixel_transpb,canalAlpha,&NotInBlackhole,&nbCollision); //Simulation du rayon principal

    int Nray=0;


    while (Nray<4 && NotInBlackhole){

        if (Nray==0){ xpixel=-(i+rdr.delta - (rdr.width + 1.0) / 2)*range;   ypixel=-(j - (rdr.height + 1.0) / 2)*range; }
        if (Nray==1){ xpixel=-(i-rdr.delta - (rdr.width + 1.0) / 2)*range;   ypixel=-(j - (rdr.height + 1.0) / 2)*range; }
        if (Nray==2){ xpixel=-(i - (rdr.width + 1.0) / 2)*range;   ypixel=-(j + rdr.delta - (rdr.height + 1.0) / 2)*range; }
        if (Nray==3){ xpixel=-(i - (rdr.width + 1.0) / 2)*range;   ypixel=-(j - rdr.delta - (rdr.height + 1.0) / 2)*range; }

        Nray++;


        xpf=xp0;
        ypf=yp0;
        zpf=zp0;

        sim_opt(x,y,z,xpixel,ypixel,&xpf,&ypf,&zpf,&NotInBlackhole); //Simulation des rayons secondaires

        if(NotInBlackhole){
            //Calcul de la taille de la zone d'impact des differents rayons (zone en terme de direction)
            //On trouve le rayon maxdist du cercle centré sur la direction du rayon principal, qui contient les autres rayons (tjrs en terme de direction)
            //Calcul en coord cartesienne pour éviter les discontinuités liées au modulo 2PI

            /*float dx0=xpcentre0-xp0; 
            float dy0=ypcentre0-yp0; 
            float dz0=zpcentre0-zp0;

            float dist2_0=sqrnorm(dx0,dy0,dz0);*/


            float dx=xpcentref-xpf;  
            float dy=ypcentref-ypf;
            float dz=zpcentref-zpf;

            float dist2=sqrnorm(dx,dy,dz);

            //if (maxdist2_0<dist2_0){ maxdist2_0=dist2_0; }
            if (maxdist2<dist2){ maxdist2=dist2; }


            
            }
        }



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


            
            //float b=dilatation*dilatation;
            float b=5e-6/maxdist2;
            *pixelr*=b;
            *pixelg*=b;
            *pixelb*=b;
        }
        

        getFinalColor(pixelr,pixelg,pixelb,pixel_transpr,pixel_transpg,pixel_transpb,canalAlpha,nbCollision);// Ajouter les superpositions liées au disque de poussière (calcul de transparence)
        
}
        


void render(const char* name){
    float* image = new float[rdr.width * rdr.height * CHANNEL_NUM]; //Tableau qui va contenir les données
    

    int Chunkcomputed=0;

    

    pl::async_par_for(0, rdr.TotalChunknumber, [&](unsigned h){
    printf("file %s - Chunk %03d started - %.2f% \n",name,h,100.*(float)Chunkcomputed/((float)rdr.TotalChunknumber));
    Chunkcomputed++;

    int j0=rdr.linesPerChunk*h;
  
     for (int j = j0; j <j0+rdr.linesPerChunk; ++j)
     {
         
        for (int i = 0; i < rdr.width; ++i)
        { 
        float x,y,z;
        float pixelr,pixelg,pixelb;   

        //Position initial du rayon à l'itération 0
        x=scn.camera.x;
        y=scn.camera.y;
        z=scn.camera.z; 

        pixelr=0.;
        pixelg=0.;
        pixelb=0.; 

        sim_bundle(i,j,x,y,z,&pixelr,&pixelg,&pixelb);

        int loc=(rdr.width*j+i)*CHANNEL_NUM;


        image[loc]   =  pixelr;
        image[loc+1] = pixelg;
        image[loc+2] = pixelb;

        } 
     }
    });


     //Postprocess possible ici

     int index = 0;
     
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
    rdr.R_inf=19.; //distance à partir de laquelle on considere etre a l'infini

    rdr.linesPerChunk=5;
    rdr.stepmax=0.009;
    rdr.stepmin=0.0005;

    rdr.maxtransparency=4;//Nombre de collision max pour un rayon



    bh.a=0.5;
    bh.precalc();//A executer avant inner orbit

    disk.R_min=bh.inner_orbit();

    disk.R_max=15.;
    disk.betamax=0.50;//Vitesse du disque au plus proche du trou noir (c'est la qu'est la vitesse max pour un profil de vitesse en /r^-1/2)
    disk.TMax=10000.;//temperature du disque au plus proche du trou noir

    scn.FOV=40.;
    
    disk.texture_rep=4.; //Répéter la texture du dique (en longeur pour ne pas qu'elle soit pixelisée)

    
    rdr.precalc(bh.a2);//quelques calculs pour avoir les carrés de certaines qtité et le fov en radian etc..
    scn.precalc(rdr.width);
    disk.precalc();

    scn.generateSky();

    //Rendre une image fixe:
    scn.camera.x=18.;
    scn.camera.y=0.;
    scn.camera.z=1.;
    scn.camera.theta=0.04;
    scn.camera.phi=0.01;

   

    printf("Start \n");
    render("resultat.png");
    printf("End \n");

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
        sphericalToCartesian(thetacoordcamera,phicoordcamera,&xcamera,&ycamera,&zcamera);
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