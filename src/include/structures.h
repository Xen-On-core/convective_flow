#ifndef BASESTRUCTURE 
#define BASESTRUCTURE

/*
 * We define the following variables describing the grid for calculations:
 *      N : number of nodes by X
 *      M : number of nodes by Y
 *      K : number of nodes by time
 *      x_0, y_0 : X and Y values at left boundaries
 *      x_1, y_1 : X and Y values at right boundaries
 *      t_1 : end time value
 *      hx, hy, tau : step by X, Y and time respectively
 *      
 * Defining constants:
 *      Re : Reynolds number
 *      Pr : Prandtl number
 *      Gr : Grashof number
 *      xi : coefficient of thermal conductivity
 *      nu : coefficient of kinematic viscosity
 *      right_beta : coefficient of thermal expansion
 */


typedef struct bounds {
    float x_0 ;
    float y_0 ;
    float x_1 ;
    float y_1 ;
    float t_1 ;
} Bounds;

typedef struct steps{
    int N;
    int M;
    int K;
    float hx;
    float hy;
    float tau;
} Steps;

typedef struct constants
{
    float Re;
    float Gr;
    float Pr;
    float xi;
    float nu;
    float rigth_beta;
} Constants;
 
typedef struct solution
{
    float *x;
    float *y;
    float *t;
    float ***temperature;
    float ***temperature12;
    float ***omega;
    float ***omega12;
    float **psi;
    float **u;
    float **v;
    float *alpha_x;
    float *beta_x;
    float *alpha_y;
    float *beta_y;
    Constants *Constants;
    Bounds *Bounds; 
    Steps *SM;
} Solution;

#endif