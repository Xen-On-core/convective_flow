#ifndef BASESTRUCTURE 
#define BASESTRUCTURE

typedef float* Vector;
typedef float** Matrix; 
typedef float*** Matrix3D; 


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
    Vector x;
    Vector y;
    Vector t;
    Matrix3D temperature;
    Matrix3D temperature12;
    Matrix3D omega;
    Matrix3D omega12;
    Matrix psi;
    Matrix u;
    Matrix v;
    Vector alpha_x;
    Vector beta_x;
    Vector alpha_y;
    Vector beta_y;
    Constants *Constants;
    Bounds *Bounds; 
    Steps *STP;
} Solution;

typedef struct starter_pack {
    int N;
    int M;
    int K;
    float x_1;
    float y_1;
    float t_1;
    float T_0;
    float T_1;
    float Re;
    float Gr;
    float Pr;
}StarterPack;

#endif