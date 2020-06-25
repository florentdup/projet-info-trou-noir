#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//Un des deux ou aucun des deux
#define DRAWSTARS 1 //Si 1 etoiles sur le fond celeste, sinon rien
#define DRAWGRID 0  //Si 1 grille sur le fond celeste

//Un des deux ou aucun des deux
#define ADISK_NORMAL 1 //Si 1, texture + temperature, sinon rien
#define ADISK_GRID 0   //Si 1, grille, sinon rien

#include "fct.h"
#include "raytracing.h"
#include "par_for.h"

#define CHANNEL_NUM 3

#define MAXITER 500000

static const int Neq = 6; //Nombre de variables pour l'integration numerique

struct Rendu rdr;
struct Scene scn;
struct Blackhole bh;
struct Disk disk;

int adisk_width, adisk_height, adisk_bpp;
uint8_t *adisk = stbi_load("adisk_upscaled.png", &adisk_width, &adisk_height, &adisk_bpp, 3);

static const int SpectrumSampleSize = 76; //nombre d'échantillons
float deltaWaveLength = 4.;               //(En nm, espacement des échantillons)
float wavelengthSamples[SpectrumSampleSize], wavelengthSamples5[SpectrumSampleSize], sensitivitySamplesR[SpectrumSampleSize], sensitivitySamplesG[SpectrumSampleSize], sensitivitySamplesB[SpectrumSampleSize];

void readSensitivityData(const char *filename, float *wavelengthSamples, float *wavelengthSamples5, float *sensitivitySamplesR, float *sensitivitySamplesG, float *sensitivitySamplesB)
{
    int cnt = 0;
    ifstream source;       // build a read-Stream
    source.open(filename); // open data

    for (std::string line; std::getline(source, line);) //read stream line by line
    {
        std::istringstream in(line); //make a stream for the line itself

        in >> wavelengthSamples[cnt] >> sensitivitySamplesR[cnt] >> sensitivitySamplesG[cnt] >> sensitivitySamplesB[cnt]; //now read the whitespace-separated floats

        wavelengthSamples5[cnt] = pow(wavelengthSamples[cnt], 5);
        cnt++;
    }
}

void getBodyColor(float *rgbR, float *rgbG, float *rgbB, float temperature, float brightness)
{
    float I;

    for (int l = 0; l < SpectrumSampleSize; ++l)
    {

        I = (6.e14 / wavelengthSamples5[l]) / exp(1.43913e7 / (wavelengthSamples[l] * temperature) - 1.) * brightness; //6e14 pour renormaliser les composante des pixels

        I *= deltaWaveLength;
        *rgbR += sensitivitySamplesR[l] * I;
        *rgbG += sensitivitySamplesG[l] * I;
        *rgbB += sensitivitySamplesB[l] * I;
    }
}

//Affiche une grille sur le disque
void getDiskColorGrid(float phi, float r, float *rgbR, float *rgbG, float *rgbB)
{
    phi = (phi + M_PI) / (2. * M_PI);

    bool a = int((100 * phi)) % 2; //100 alternances de couleur par tour

    bool b = ((r - disk.R_min) / (disk.R_max - disk.R_min) > .5); //Séparer le disque en 2 radialement

    if (a ^ b) //Ou exclusif
    {
        *rgbR = 255.;
        *rgbG = 0.;
    }
    else
    {
        *rgbR = 0.;
        *rgbG = 255.;
    }
}

//Passage aux coordonnees de Boyer Lindquist
//Pour l'instant juste passage au coordonnées sphériques (plus simple)
//Utilisé uniquement pour trouver la pos initiale de la camera, pas important
void cartesianToBl(float x, float y, float z, float *r, float *theta, float *phi)
{

    float r2 = x * x + y * y + z * z;
    *r = sqrt(r2);

    *phi = atan2(y, x);
    *theta = acos(z / (*r));
}

//Passage des coordonnées de Boyer Lindquist au coordonnées cartésienne (exacte)
void blToCartesian(float r, float theta, float phi, float *x, float *y, float *z)
{

    float sintheta = sin(theta);
    float costheta = cos(theta);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float temp = sintheta * sqrt(r * r + bh.a2);

    *x = temp * cosphi;
    *y = temp * sinphi;
    *z = r * costheta;
}

/* Fonction utilisé pour l'intégration */
void geodesic(float L, float kappa, float *y, float *dydx)
{
    float r, theta, pr, ptheta;

    r = y[0];
    theta = y[1];
    pr = y[4];
    ptheta = y[5];

    float r2 = r * r;
    float twor = 2.0 * r;

    float sintheta, costheta;
    sintheta = sin(theta);
    costheta = cos(theta);

    float costheta2 = costheta * costheta;
    float sintheta2 = sintheta * sintheta;

    float sigma = r2 + bh.a2 * costheta2;
    float delta = r2 - twor + bh.a2;
    float sd = sigma * delta;
    float siginv = 1.0 / sigma;
    float bot = 1.0 / sd;

    /* Prevent problems with the axis */
    if (sintheta < 1e-8)
    {
        sintheta = 1e-8;
        sintheta2 = 1e-16;
    }

    dydx[0] = -pr * delta * siginv;
    dydx[1] = -ptheta * siginv;
    dydx[2] = -(twor * bh.a + (sigma - twor) * L / sintheta2) * bot;
    dydx[3] = -(1.0 + (twor * (r2 + bh.a2) - twor * bh.a * L) * bot);
    dydx[4] = -(((r - 1.0) * (-kappa) + twor * (r2 + bh.a2) - 2.0 * bh.a * L) * bot - 2.0 * pr * pr * (r - 1.0) * siginv);
    dydx[5] = -sintheta * costheta * (L * L / (sintheta2 * sintheta2) - bh.a2) * siginv;
}

/* Conditions initiales pour un rayon */
void initial(float r0, float theta0, float *L, float *kappa, float *y0, float *ydot0, float x, float y)
{
    y0[0] = r0;
    y0[1] = theta0;
    y0[2] = 0;
    y0[3] = 0;
    y0[4] = cos(y) * cos(x);
    y0[5] = sin(y) / r0;

    float sintheta, costheta;
    sintheta = sin(theta0);
    costheta = cos(theta0);
    float costheta2 = costheta * costheta;
    float sintheta2 = sintheta * sintheta;

    float rdot0 = y0[4];
    float thetadot0 = y0[5];

    float r2 = r0 * r0;
    float sigma = r2 + bh.a2 * costheta2;
    float delta = r2 - 2.0 * r0 + bh.a2;
    float s1 = sigma - 2.0 * r0;

    y0[4] = rdot0 * sigma / delta;
    y0[5] = thetadot0 * sigma;

    ydot0[0] = rdot0;
    ydot0[1] = thetadot0;
    ydot0[2] = cos(y) * sin(x) / (r0 * sin(theta0));

    float phidot0 = ydot0[2];
    float energy2 = s1 * (rdot0 * rdot0 / delta + thetadot0 * thetadot0) + delta * sintheta2 * phidot0 * phidot0;

    float energy = sqrt(energy2);

    /* Energie de 1 */
    y0[4] = y0[4] / energy;
    y0[5] = y0[5] / energy;

    /* Angular Momentum with E = 1 */
    *L = ((sigma * delta * phidot0 - 2.0 * bh.a * r0 * energy) * sintheta2 / s1) / energy;

    *kappa = y0[5] * y0[5] + bh.a2 * sintheta2 + (*L) * (*L) / sintheta2;

    /* Hack - make sure everything is normalized correctly by a call to geodesic */

    geodesic(*L, *kappa, y0, ydot0);
}

//Obtenir la direction du rayon à partir des coord de BL (y) et de leurs dérivées (dydx)
//Retourne un vecteur normé
void getDirection(float *y, float *dydx, float *xp, float *yp, float *zp)
{

    float costheta = cos(y[1]);
    float sintheta = sin(y[1]);
    float cosphi = cos(y[2]);
    float sinphi = sin(y[2]);
    float r2 = y[0] * y[0];
    float R2 = r2 + bh.a2;
    float R = sqrt(r2 + bh.a2);

    *xp = R * dydx[1] * cosphi * costheta - R * dydx[2] * sinphi * sintheta + sintheta * cosphi * y[0] * dydx[0] / R;
    *yp = R * dydx[1] * sinphi * costheta + R * dydx[2] * cosphi * sintheta + sintheta * sinphi * y[0] * dydx[0] / R;
    *zp = -y[0] * sintheta * dydx[1] + dydx[0] * costheta;
    normalise(xp, yp, zp);
}

//Simulation complete d'un rayon (trajectoire+collision)
void sim(float r0, float theta0, float phi0, float xpixel, float ypixel, float *xp, float *yp, float *zp, float *pixel_transpr, float *pixel_transpg, float *pixel_transpb, float *canalAlpha, bool *ReachedInfinity, int *nbCollision)
{
    int N = 0; //Nombre d'itérations
    int k = 0; //nombre de collision

    float oldtheta; //pour detecter le passage à travers le plan z=0 (pour tracer le disque)

    bool zSignChange, diskDistance, diskCollision;

    //Euler

    //float y[Neq];
    //float dydx[Neq];

    //RK4
    float y[Neq];
    float ak[Neq];
    float dydx1[Neq], dydx2[Neq], dydx3[Neq], dydx4[Neq];
    float ytemp[Neq];

    float L, kappa;

    initial(r0, theta0, &L, &kappa, y, ak, xpixel, ypixel);

    float currentStep = rdr.step(y[0]);

    int l = 0; //Juste pour effectuer des boucles sans redeclarer de var

    while ((N < MAXITER) && (rdr.R_min < y[0]) && (y[0] < rdr.R_inf))
    {

        N += 1;

        oldtheta = y[1];

        currentStep = rdr.step(y[0]);

        //Euler
        /*for (int l = 0; l < Neq; l++)
	    {
		    float hdydx = currentStep * dydx[l];
		    y[l] = y[l] + hdydx;
	    }

	    geodesic(L,kappa,y, dydx);*/

        //RK4
        geodesic(L, kappa, y, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx1[l] = ak[l];
            ytemp[l] = y[l] + .5 * currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx2[l] = ak[l];
            ytemp[l] = y[l] + .5 * currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx3[l] = ak[l];
            ytemp[l] = y[l] + currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx4[l] = ak[l];
            dydx1[l] = currentStep / 6. * (dydx1[l] + 2. * dydx2[l] + 2. * dydx3[l] + dydx4[l]); //Recuperer une bonne approx de la derivée dans dydx1
            y[l] = y[l] + dydx1[l];
        }

        zSignChange = (oldtheta > M_PI_2) != (y[1] > M_PI_2);      //on traverse le plan z=0
        diskDistance = (y[0] < disk.R_max) && (y[0] > disk.R_min); //On l'a traversé la ou est le disque
        diskCollision = zSignChange && diskDistance;

        if (diskCollision)
        {
            float xpos, ypos, zpos;
            blToCartesian(y[0], y[1], y[2], &xpos, &ypos, &zpos); //Obtenir la position en coord cartesiennes
            getDirection(y, dydx1, xp, yp, zp);                   //Obtenir la direction (en coord cartesiennes)

            //Trouver le point de colision (en allant tout droit entre les deux point au dessus et en dessous du disque (valable si le pas est assez petit))
            float lambda = -zpos / (*zp);
            float coll_x = xpos + lambda * (*xp);
            float coll_y = ypos + lambda * (*yp);
            float coll_z = zpos + lambda * (*zp);
            float r = norm(coll_x, coll_y, coll_z);

            diskCollision = (r < disk.R_max) && (r > disk.R_min); //#reverification plus précise
            if (diskCollision)
            {
                if (k < rdr.maxtransparency)
                {

                    float phi = atan2(coll_y / r, coll_x / r); //Coordonné du point d'impact (en sphérique)

#if ADISK_NORMAL == 1

                    float diskspeed_x, diskspeed_y, diskspeed_z;

                    //rotation prograde du disque

                    sphericalToCartesian(M_PI / 2., phi + M_PI / 2., &diskspeed_x, &diskspeed_y, &diskspeed_z); //Direction de la vitesse en ce point, sens trigo,mvt circulaire

                    float beta = disk.RotationSpeed(r);                                                         //Obtenir la norme de la vitesse des poussières en ce point
                    float costhetadoppler = -((*xp) * diskspeed_x + (*yp) * diskspeed_y + (*zp) * diskspeed_z); //Obtenir le cosinus de l'angle entre le rayon et la vitesse des particules (Les vecteurs sont normés)
                    //- à cause du sens des rayons

                    //Temperature du disque à cet endroit
                    float Temp = disk.Temp(r);

                    //Décalage en fréquence égal à décalage en temperature, valable aussi pour l'intensité

                    Temp = Temp * (1 + beta * costhetadoppler) / sqrt(1 - beta * beta); //Effet doppler relativiste

                    //calcul du redshift gravitationnel
                    float gtt = 1. - 1. / r;
                    float gtphi = bh.a / r;
                    float gphiphi = -(r * r + bh.a2 + bh.a2 / r);
                    float omega = beta / r;                                                //Vitesse angulaire
                    Temp = Temp * sqrt(gtt + 2 * gtphi * omega + gphiphi * omega * omega); //Redshift gravitationnel

                    //Recuperer la texture du disque à cet endroit, qu'on utilise uniquemenet pour obtenir la transparence du disque
                    int cx = int((r - disk.R_min) * (adisk_height - 1) / (disk.R_max - disk.R_min));
                    int cy = int((adisk_width - 1) * mod(disk.texture_rep * phi - M_PI, 2. * M_PI) / (2. * M_PI));

                    int loc = (cx * adisk_width + cy) * CHANNEL_NUM;

                    pixel_transpr[k] = 0.;
                    pixel_transpg[k] = 0.;
                    pixel_transpb[k] = 0.;

                    getBodyColor(&pixel_transpr[k], &pixel_transpg[k], &pixel_transpb[k], Temp, 1.); //Obtenir la couleur

                    canalAlpha[k] = (float)adisk[loc + 1] / 255.; //Choix: on utilise la composante verte (le +1) pour reconstituer la transparence

#endif

#if ADISK_GRID == 1

                    getDiskColorGrid(phi, r, &pixel_transpr[k], &pixel_transpg[k], &pixel_transpb[k]); //Motif de grille sur le disque et calcul de la transparence
                    canalAlpha[0] = 1.;
#endif

                    k++;
                }
            }
        }
    }

    *ReachedInfinity = (y[0] >= rdr.R_inf);
    *nbCollision = k;

    if (*ReachedInfinity)
    {
        //Pour visualiser quels rayons ont atteints R_inf
        //*nbCollision=1;
        //pixel_transpr[0]=255;
        //canalAlpha[0]=.5;

        getDirection(y, dydx1, xp, yp, zp);
        normalise(xp, yp, zp); //Ne pas effectuer inutilement cette opération
    }
}

//Simulation d'un rayon (trajectoire uniquement)
void sim_opt(float r0, float theta0, float phi0, float xpixel, float ypixel, float *xp, float *yp, float *zp, bool *ReachedInfinity)
{

    int N = 0;

    //Euler
    //float y[Neq];
    //float dydx[Neq];

    //RK4
    float y[Neq], ak[Neq];
    float dydx1[Neq], dydx2[Neq], dydx3[Neq], dydx4[Neq];
    float ytemp[Neq];

    float L, kappa;

    initial(r0, theta0, &L, &kappa, y, ak, xpixel, ypixel); //CI

    float currentStep = rdr.step(y[0]);

    int l;

    while ((N < MAXITER) && (rdr.R_min < y[0]) && (y[0] < rdr.R_inf))
    {

        N += 1;

        currentStep = rdr.step(y[0]);

        //Euler
        /*for (int l = 0; l < Neq; l++)
	    {
		    float hdydx = currentStep * dydx[l];
		    y[l] = y[l] + hdydx;
	    }
	    geodesic(L,kappa,y, dydx);*/

        //RK4

        geodesic(L, kappa, y, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx1[l] = ak[l];
            ytemp[l] = y[l] + .5 * currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx2[l] = ak[l];
            ytemp[l] = y[l] + .5 * currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx3[l] = ak[l];
            ytemp[l] = y[l] + currentStep * ak[l];
        }
        geodesic(L, kappa, ytemp, ak);
        for (l = 0; l < Neq; ++l)
        {
            dydx4[l] = ak[l];
            dydx1[l] = currentStep / 6. * (dydx1[l] + 2. * dydx2[l] + 2. * dydx3[l] + dydx4[l]); //Recuperer une bonne aprrox de la derivée
            y[l] = y[l] + dydx1[l];
        }
    }

    *ReachedInfinity = (y[0] >= rdr.R_inf);

    if (*ReachedInfinity)
    { //Ne pas effectuer inutilement cette opération
        getDirection(y, dydx1, xp, yp, zp);
        normalise(xp, yp, zp);
    }
}

//Simulation d'un faisceu de rayons (uniquement dans le cas ou on dessine des étoiles)
void sim_bundle(int i, int j, float x, float y, float z, float *pixelr, float *pixelg, float *pixelb)
{ //arguments: i et j position du pixel, x,y,z position initiale du rayon dans l'espace

    float r0, theta0, phi0;
    cartesianToBl(x, y, z, &r0, &theta0, &phi0); //Approximativement (compliqué de passer de cartesien à BL), mais pas important

    float xpcentre0, ypcentre0, zpcentre0;
    float xpcentref, ypcentref, zpcentref;
    float xp0, yp0, zp0;
    float xpf, ypf, zpf;

    bool ReachedInfinity = true;
    int nbCollision;

    float maxdist2 = 0.;
    float maxdist2_0 = 0.;
    //float mindist2=1000.;
    //float mindist2_0=1000.;

    /*float Distances2[4];
    float Distances2_0[4];
    float dilatation=1.;*/

    float pixel_transpr[rdr.maxtransparency], pixel_transpg[rdr.maxtransparency], pixel_transpb[rdr.maxtransparency];
    float canalAlpha[rdr.maxtransparency];

    xpcentref = xpcentre0;
    ypcentref = ypcentre0;
    zpcentref = zpcentre0;

    float range = 0.075 * 20. / (rdr.width - 1.0); //FOV
    //float range = 0.1 * 20. / (rdr.width - 1.0);

    float xpixel = -(i - (rdr.width + 1.0) / 2) * range;
    float ypixel = -(j - (rdr.height + 1.0) / 2) * range;

    sim(r0, theta0, phi0, xpixel, ypixel, &xpcentref, &ypcentref, &zpcentref, pixel_transpr, pixel_transpg, pixel_transpb, canalAlpha, &ReachedInfinity, &nbCollision); //Simulation du rayon principal

#if DRAWSTARS == 1

    int Nray = 0;

    while (Nray < 4 && ReachedInfinity)
    {

        if (Nray == 0)
        {
            xpixel = -(i + rdr.delta - (rdr.width + 1.0) / 2) * range;
            ypixel = -(j - (rdr.height + 1.0) / 2) * range;
        }
        if (Nray == 1)
        {
            xpixel = -(i - rdr.delta - (rdr.width + 1.0) / 2) * range;
            ypixel = -(j - (rdr.height + 1.0) / 2) * range;
        }
        if (Nray == 2)
        {
            xpixel = -(i - (rdr.width + 1.0) / 2) * range;
            ypixel = -(j + rdr.delta - (rdr.height + 1.0) / 2) * range;
        }
        if (Nray == 3)
        {
            xpixel = -(i - (rdr.width + 1.0) / 2) * range;
            ypixel = -(j - rdr.delta - (rdr.height + 1.0) / 2) * range;
        }

        Nray++;

        xpf = xp0;
        ypf = yp0;
        zpf = zp0;

        sim_opt(r0, theta0, phi0, xpixel, ypixel, &xpf, &ypf, &zpf, &ReachedInfinity); //Simulation des rayons secondaires

        if (ReachedInfinity)
        {
            //Calcule de la taille de la zone d'impact des differents rayons (zone en terme de direction)
            //On trouve le rayon maxdist du cercle centré sur la direction du rayon principal, qui contient les autres rayons (tjrs en terme de direction)
            //Calcul en coord cartesienne pour éviter les discontinuités liées au modulo 2PI

            //float dx0=xpcentre0-xp0;
            //float dy0=ypcentre0-yp0;
            //float dz0=zpcentre0-zp0;

            //float dist2_0=sqrnorm(dx0,dy0,dz0);

            float dx = xpcentref - xpf;
            float dy = ypcentref - ypf;
            float dz = zpcentref - zpf;

            float dist2 = sqrnorm(dx, dy, dz);

            //if (maxdist2_0<dist2_0){ maxdist2_0=dist2_0; }
            if (maxdist2 < dist2)
            {
                maxdist2 = dist2;
            }
        }
    }

    if (ReachedInfinity)
    { //Si tous les rayons se sont échappés, on regarde quels étoiles sont dans le cercle formé par le faisceau sur le ciel

        float maxdist = sqrt(maxdist2);
        float sqrdistanceToStar; //sqrnorm(scn.starx[i]-xpcentref,scn.stary[i]-ypcentref,scn.starz[i]-zpcentref);

        float c1, c2, c3;

        for (int i = 0; i < Nstar; ++i)
        {

            c1 = abs(scn.starx[i] - xpcentref);

            if (c1 < maxdist)
            {
                c2 = abs(scn.stary[i] - ypcentref);
                if (c2 < maxdist)
                {
                    c3 = abs(scn.starz[i] - zpcentref);
                    if (c3 < maxdist)
                    {
                        sqrdistanceToStar = sqrnorm(c1, c2, c3);
                        if (sqrdistanceToStar < maxdist2)
                        {
                            getBodyColor(pixelr, pixelg, pixelb, scn.starTemp[i], scn.starBrightness[i] * 12.); //(Récupere la couleur et effectue l'addition composante par composante)
                        }
                    }
                }
            }
        }

        float b = 1e-7 / maxdist2; //Angle solide initial environ constant égal à 1e-7
        *pixelr *= b;
        *pixelg *= b;
        *pixelb *= b;
    }
#endif
#if DRAWGRID == 0
    getFinalColor(pixelr, pixelg, pixelb, pixel_transpr, pixel_transpg, pixel_transpb, canalAlpha, nbCollision); // Calcul de transparence
#else
    getFinalColorGrid(xpcentref, ypcentref, zpcentref, pixelr, pixelg, pixelb, pixel_transpr, pixel_transpg, pixel_transpb, nbCollision, ReachedInfinity); // Calcul de la transparence, et calcul de la couleur de la grille
#endif
}

void render(const char *name)
{
    float *image = new float[rdr.width * rdr.height * CHANNEL_NUM]; //Tableau qui va contenir les données

    int Chunkcomputed = 0;

    pl::async_par_for(0, rdr.TotalChunknumber, [&](unsigned h) {
        //for (int h=0;h<rdr.TotalChunknumber;++h){ //A utiliser si jamais pl async ne marche pas
        cout<<"file "<< name <<"- Chunk "<< h <<"/"<< rdr.TotalChunknumber<<" started - "<< 100. * (float)Chunkcomputed / ((float)rdr.TotalChunknumber) <<"%"<<endl;  
        Chunkcomputed++;

        int j0 = rdr.linesPerChunk * h;

        for (int j = j0; j < j0 + rdr.linesPerChunk; ++j)
        {

            for (int i = 0; i < rdr.width; ++i)
            {

                //if ((i>.3*rdr.width) && (i<0.55*rdr.width) && (j>0.1*rdr.height) && (j<0.55*rdr.height)){//calculer des petites portions pour faire des tests
                //if ((i>.0*rdr.width) && (i<0.05*rdr.width) && (j>0.*rdr.height) && (j<0.05*rdr.height)){

                float pixelr, pixelg, pixelb;

                pixelr = 0.;
                pixelg = 0.;
                pixelb = 0.;

                sim_bundle(i, j, scn.camera.x, scn.camera.y, scn.camera.z, &pixelr, &pixelg, &pixelb);

                int loc = (rdr.width * j + i) * CHANNEL_NUM;

                image[loc] = pixelr;
                image[loc + 1] = pixelg;
                image[loc + 2] = pixelb;

                //}
            }
        }
    });
    //}

    //Postprocess possible ici

    int index = 0;

    uint8_t *image_byte = new uint8_t[rdr.width * rdr.height * CHANNEL_NUM]; //Image finale

    for (int j = 0; j < rdr.height; ++j)
    {
        for (int i = 0; i < rdr.width; ++i)
        {

            image_byte[index] = char(image[index]);
            index++;

            image_byte[index] = char(image[index]);
            index++;

            image_byte[index] = char(image[index]);
            index++;
        }
    }

    stbi_write_png(name, rdr.width, rdr.height, CHANNEL_NUM, image_byte, rdr.width * CHANNEL_NUM);
}

int image()
{   
    /*ifstream monFlux(fichier) ;
    if (monFlux)
    {  
        monFlux>>bh.a>>RAdisk?>>scn.camera.x>>scn.camera.y>>scn.camera.z>>rdr.R_inf>>disk.R_max>>rdr.height>>rdr.width; */
    if (adisk==NULL)
    {
    exit(0);
    }
    rdr.height = 108;
    rdr.width = 192;
    rdr.R_inf = 21.5; //distance à partir de laquelle on considere etre a l'infini

    rdr.linesPerChunk = 5; //Blocs de 5 ligne traités en parallele

    //Euler
    /*rdr.stepmax=0.0035;
    rdr.stepmin=0.002; //0.001 min sinon erreurs arrondi ?*/

    //RK4
    rdr.stepmax = 0.02;
    rdr.stepmin = 0.007;

    rdr.delta = .5; //Ecart angulaire (en pixel) entre les rayons d'un même faisceau

    bh.a = 0.5;   //spin (adimensionné,entre 0 et 1)
    bh.precalc(); //A executer avant inner orbit

    disk.R_min = bh.inner_orbit(); //Trouver l'orbite la plus proche du trou noir encore stable (prograde)

    disk.R_max = 16.;
    disk.betamax = 0.65; //Vitesse du disque au plus proche du trou noir (c'est la qu'est la vitesse max pour un profil de vitesse en r^-1/2)
    disk.TMax = 17500.;  //temperature du disque au plus proche du trou noir

    //scn.FOV=40.;

    disk.texture_rep = 8.; //Répéter la texture du dique (en longeur pour ne pas qu'elle soit pixelisée)

    rdr.precalc(bh.a2); //quelques calculs pour avoir les carrés de certaines qtité et le fov en radian etc..
    //scn.precalc(rdr.width);
    disk.precalc();

    scn.generateSky();

    //Rendre une image fixe:
    scn.camera.x = 20.;
    scn.camera.y = 0.;
    scn.camera.z = 1.2;

    /*scn.camera.theta=0.04;
    scn.camera.phi=0.01;*/

    readSensitivityData("sensitivity.txt", wavelengthSamples, wavelengthSamples5, sensitivitySamplesR, sensitivitySamplesG, sensitivitySamplesB);

    /*for (int l = 0; l < SpectrumSampleSize; ++l)
    {
        printf("WL:%f R:%f, G:%f B:%f \n",wavelengthSamples[l],sensitivitySamplesR[l],sensitivitySamplesG[l],sensitivitySamplesB[l]);
    }*/

    //float R,G,B;
    //getBodyColor(&R,&G,&B,6000,1.);
    //printf("R:%f G:%f, B:%f \n",R,G,B);

    cout<<"Start"<<endl;
    render("resultat.png");
    cout<<"End"<<endl;

    stbi_image_free(adisk);

    
    /*}
    else
    {
        cout<<"Erreur le fichier ne s'est pas ouvert"<<endl;
    }*/
    return 0;

}

/*int main()
{
    image();
    return 0;
}*/

#include<boost/python.hpp> 

BOOST_PYTHON_MODULE(kerr)
{
    using namespace boost::python;
    def("image", image);
}