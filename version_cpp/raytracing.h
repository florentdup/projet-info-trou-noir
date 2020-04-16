#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>

struct Camera {
	float x, y, z,theta,phi;
};


struct Scene {
	Camera camera;
    float FOV;
    float RAdisk_min;
    float RAdisk_max;
    float adisk_texture_rep;

    float RAdisk_min2;
    float RAdisk_max2;

    float fov_rad;
    float distance_viewplane;


    void initScene(int width){
        RAdisk_min2=pow((RAdisk_min-0.5),2); //#le -1 et +1 sont une marge pour la detection avant de calculer le point de collision exact
        RAdisk_max2=pow((RAdisk_max+0.5),2)+2.;
        fov_rad=FOV*M_PI/180.;
        distance_viewplane=width/(2*tan(fov_rad/2.));
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
    
    
    void initRendering(){
        R_min=R_schwarzschild;
        R_inf2=pow(R_inf,2);
        R_min2=pow(R_min,2);
    }


};
