#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>

float FloatRand(float MaxVal)
{
    return ((float)rand() / (float)RAND_MAX) * MaxVal;
}

void cartesianToSpherical(float x, float y, float z, float *theta, float *phi)
{
    *theta = acos(z);
    *phi = atan2(y, x);
}

void sphericalToCartesian(float theta, float phi, float *x, float *y, float *z)
{
    float tmp = sin(theta);
    *x = tmp * cos(phi);
    *y = tmp * sin(phi);
    *z = cos(theta);
}

float mod(float x, float y)
{
    int resultat = floor(x / y);
    return x - (resultat * y);
}

float sqrnorm(float x, float y, float z)
{
    return x * x + y * y + z * z;
}
float norm(float x, float y, float z)
{
    return sqrt(x * x + y * y + z * z);
}

void normalise(float *x, float *y, float *z)
{
    float tmp = sqrt(sqrnorm(*x, *y, *z));
    *x /= tmp;
    *y /= tmp;
    *z /= tmp;
}

float avg(float *array, int size)
{
    float res = 0;
    for (int i = 0; i < size; ++i)
    {
        res += array[i];
    }
    res = res / (float)size;
    return res;
}

void normalisepixel(float *rgbR, float *rgbG, float *rgbB)
{
    float maxcomp = *rgbR;
    if (*rgbG > maxcomp)
    {
        maxcomp = *rgbG;
    }
    if (*rgbB > maxcomp)
    {
        maxcomp = *rgbB;
    }

    if (maxcomp > 255.)
    {
        *rgbR = 255. * (*rgbR) / maxcomp;
        *rgbG = 255. * (*rgbG) / maxcomp;
        *rgbB = 255. * (*rgbB) / maxcomp;
    }
}

//Utilisé pour ajouter à l'image du fond celest le disque de poussière en superposant dans l'ordre les masques
void getFinalColor(float *rgbR, float *rgbG, float *rgbB, float *pixel_transpr, float *pixel_transpg, float *pixel_transpb, float *canalAlpha, int nbCollision)
{

    normalisepixel(rgbR, rgbG, rgbB);

    for (int l = nbCollision - 1; l >= 0; --l)
    {

        normalisepixel(&pixel_transpr[l], &pixel_transpg[l], &pixel_transpb[l]);

        *rgbR = canalAlpha[l] * pixel_transpr[l] + (1 - canalAlpha[l]) * (*rgbR);
        *rgbG = canalAlpha[l] * pixel_transpg[l] + (1 - canalAlpha[l]) * (*rgbG);
        *rgbB = canalAlpha[l] * pixel_transpb[l] + (1 - canalAlpha[l]) * (*rgbB);
    }
}

void getFinalColorGrid(float xpcentref, float ypcentref, float zpcentref, float *rgbR, float *rgbG, float *rgbB, float *pixel_transpr, float *pixel_transpg, float *pixel_transpb, int nbCollision, bool ReachedInfinity)
{
    float theta, phi;
    cartesianToSpherical(xpcentref, ypcentref, zpcentref, &theta, &phi);
    phi = (phi + M_PI) / (2. * M_PI);
    theta /= M_PI;

    bool a = int((100 * theta)) % 2; //100 changements par tour
    bool b = int((100 * phi)) % 2;   //100 changements par tour

#if ADISK_GRID == 1

    if (nbCollision > 0)
    {
        *rgbR = pixel_transpr[0];
        *rgbG = pixel_transpg[0];
        *rgbB = pixel_transpb[0];
    }
    else if (ReachedInfinity)
    {
        *rgbR = b * 200;
        *rgbG = 0.;
        *rgbB = a * 200;
    }
#else

    if (ReachedInfinity)
    {
        *rgbR = b * 200;
        *rgbG = 0.;
        *rgbB = a * 200;
    }
#endif
}
