#ifndef UTILS_H
#define UTILS_H

/*
 * Defining constants:
 *      Re : Reynolds number
 *      Pr : Prandtl number
 *      Gr : Grashof number
 *      xi : coefficient of thermal conductivity
 *      nu : coefficient of kinematic viscosity
 *      right_beta : coefficient of thermal expansion
 */
typedef struct ConvectiveFlowParameters {
    double Re;
    double Gr;
    double Pr;
    double xi;
    double nu;
    double rigth_beta;
    double lambda;
} CFP;


#endif /* UTILS_H */