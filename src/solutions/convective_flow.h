#ifndef SOLUTIONS_H
#define SOLUTIONS_H

#include "array.h"
#include "grid.h"

typedef struct ConvectiveFlow1D
{
    Grid1D grid;
    Vector temperature;
    Vector omega;
    Vector psi;
    Vector u;
} CFS_1D;

typedef struct ConvectiveFlow2D
{
    Grid2D grid;
    Matrix temperature;
    Matrix omega;
    Matrix psi;
    Matrix u;
    Matrix v;
} CFS_2D;

typedef struct ConvectiveFlow3D
{
    Grid3D grid;
    Matrix3D temperature;
    Matrix3D omega;
    Matrix3D psi;
    Matrix3D u;
    Matrix3D v;
    Matrix3D w;
} CFS_3D;

void convective_flow_solution_1d();

void convective_flow_solution_2d();

void convective_flow_solution_3d();

#endif /* SOLUTIONS_H */