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

    float step;

    float R_schwarzschild;

    float R_inf;
    float R_min;
    
    float R_inf2;
    float R_min2;

    int maxtransparency;

    
    void precalc(){
        R_min=R_schwarzschild;
        R_inf2=R_inf*R_inf;
        R_min2=R_min*R_min;
    }


};


struct blackhole {
    float a;
    float R;
    float M;
    
};


struct disk {
    
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
