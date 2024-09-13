#ifndef GRID_H
#define GRID_H

typedef enum BoundaryConditionType
{
    DIRICHLET,
    NEUMANN,
    ROBIN,
    CAUCHY
} BCType;

#endif /* GRID_H */