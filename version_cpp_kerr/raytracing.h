#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>




const int Nstar=50000;//Nombre d'étoiles

struct Camera {
	float x, y, z,theta,phi;
};


struct Scene {

    
	Camera camera;
    float FOV;
    
    

    float fov_rad;
    float distance_viewplane;

    float starx[Nstar],stary[Nstar],starz[Nstar],starTemp[Nstar],starBrightness[Nstar];

    
    void precalc(int width){
        fov_rad=FOV*M_PI/180.;
        distance_viewplane=width/(2*tan(fov_rad/2.));
    }
    

    void generateSky(){

    for (int i = 0; i<Nstar; ++i)
    {
        float thetastar=FloatRand(M_PI);
        float phistar=FloatRand(2*M_PI)-M_PI;

        starTemp[i]=3000+FloatRand(5000); //Temperature de l'etoile en Kelvin
        starBrightness[i]=.9+FloatRand(.1); //Fluctuations de luminosité des étoiles
        sphericalToCartesian(thetastar,phistar,&starx[i],&stary[i],&starz[i]); //Coordonnées des etoiles sur la sphere unitée cartesienne

    }
}
    


};

struct Rendu {
    int height;
    int width;

    float defaultstep;
    float R_schwarzschild;

    float R_inf;
    float R_min;
    
    float R_inf2;
    float R_min2;

    int maxtransparency;

    int linesPerChunk,TotalChunknumber;

    

    float C;

    float delta=.9;//"Taille angulaire (en pixel)" sur la direction des rayons d'un même groupe  

    float stepmax;
    float stepmin;
    float dstep;

    float step(float r){
        float res=stepmin+dstep*C*(r-R_min)*(r-R_min);
        return res;
    }

    
    void precalc(float a2){
        
        R_min=1.0 + sqrt(1.0-a2) + 1e-5;//Horizon du trou noir

        R_inf2=R_inf*R_inf;
        R_min2=R_min*R_min;

        TotalChunknumber=height/linesPerChunk;

        dstep=stepmax-stepmin;
        C=1./((R_inf-R_min)*(R_inf-R_min));
    }


};


struct Disk {
    
    float TMax;

    float R_min, R_max;
    float texture_rep;

    float betamax;//Vitesse de roation maximale du disque (au centre)

    float R_min2;
    float R_max2;

    float A,B;//Ceofficient pour le calcul de la temperature du disque

    float RotationSpeed(float r){
        return B/sqrt(r);  //Rotation en r^-1/2
    }

    float Temp(float r){
        return A*pow(r,-0.75) ;
    }
    void precalc(){
        R_min2=pow((R_min-0.5),2); //#le -1 et +1 sont une marge pour la detection avant de calculer le point de collision exact
        R_max2=pow((R_max+0.5),2)+2.;
        A=TMax/pow(R_min,-0.75);  //Defini pour fixer la temperature de couleur atteinte en point le plus chaud du disque (le plus proche du trou noir)
        B=sqrt(R_min)*betamax;
    }

    
};



struct Blackhole {
    float a,a2;

    float inner_orbit(void){
	    float z1 = 1+cbrt(1-a2)*(cbrt(1+a)+cbrt(1-a));//racine troisième
	    float z2 = sqrt(3*a2+z1*z1);
	    return 3+z2-sqrt((3-z1)*(3+z1+2*z2));
    }
    

    void precalc(){
        a2*a*a;
    }
    
};




