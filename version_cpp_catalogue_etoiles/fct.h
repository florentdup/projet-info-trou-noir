#include <stdint.h>
#include <iostream>
#include <math.h>
#include <tgmath.h>



float FloatRand( float MaxVal )
{
    return ( (float)rand( ) / (float)RAND_MAX ) * MaxVal;
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

float avg(float* array,int size){
    float res=0;
    for (int i=0;i<size;++i){
        res+=array[i];
    }
    res=res/(float)size;
    return res;
}

