#ifndef GRID_H
#define GRID_H

#include "array.h"


typedef struct Grid1D
{
    Vector *x1;
    float x1_l;
    float x1_r;
    double hx1;
} Grid1D;

typedef struct Grid2D
{
    Vector *x1;
    Vector *x2;
    float x1_l;
    float x2_l;
    float x1_r;
    float x2_r;
    double hx1;
    double hx2;
} Grid2D;

typedef struct Grid3D
{
    Vector *x1;
    Vector *x2;
    Vector *x3;
    float x1_l;
    float x2_l;
    float x3_l;
    float x1_r;
    float x2_r;
    float x3_r;
    double hx1;
    double hx2;
    double hx3;
} Grid3D;

typedef struct PolarGrid
{
    Vector rho;
    Vector phi;
    float rho_l;
    float phi_l;
    float rho_r;
    float phi_r;
    double hrho;
    double hphi;
} PolarGrid;

typedef struct CylindricalGrid
{
    Vector rho;
    Vector phi;
    Vector z;
    float rho_l;
    float phi_l;
    float z_l;
    float rho_r;
    float phi_r;
    float z_r;
    double hrho;
    double hphi;
    double hz;
} CylindricaGrid;

typedef struct SphericalGrid
{
    Vector rho;
    Vector phi;
    Vector theta;
    float rho_l;
    float phi_l;
    float theta_l;
    float rho_r;
    float phi_r;
    float theta_r;
    double hrho;
    double hphi;
    double htheta;
} SphericalGrid;

void get_grid_info();

#endif /* GRID_H */
