#ifndef GRID_H
#define GRID_H

#include "array.h"


typedef struct Grid1D
{
    int size;
    double hx1;
    Vector x1;
} Grid1D;

typedef struct Grid2D
{
    int x1_size;
    int x2_size;
    float x1_r;
    float x2_r;
    float x1_l;
    float x2_l;
    double hx1;
    double hx2;
    Vector x1;
    Vector x2;
} Grid2D;

typedef struct Grid3D
{
    int x1_size;
    int x2_size;
    int x3_size;
    float x1_r;
    float x2_r;
    float x3_r;
    float x1_l;
    float x2_l;
    float x3_l;
    double hx1;
    double hx2;
    double hx3;
    Vector x1;
    Vector x2;
    Vector x3;
} Grid3D;

typedef struct PolarGrid
{
    int rho_size;
    int phi_size;
    Vector rho;
    Vector phi;
} PolarGrid;

typedef struct CylindricalGrid
{
    int rho_size;
    int phi_size;
    int z_size;
    Vector rho;
    Vector phi;
    Vector z;
} CylindricaGrid;

typedef struct SphericalGrid
{
    int rho_size;
    int phi_size;
    int theta_size;
    int psi_size;
    Vector rho;
    Vector phi;
    Vector theta;
    Vector psi;
} SphericalGrid;



#endif /* GRID_H */