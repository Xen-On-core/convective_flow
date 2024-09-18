#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include "utils/array.h"

float epsilon = 1e-6;
#define NUMCHECK(number) \
            if (abs(number - 0) < epsilon) { \
                printf("Wrong value of %s: %f\n", #number, number); \
                exit(1); \
            }

/*
 * We define the following variables describing the grid for calculations:
 *      N : number of nodes by X
 *      M : number of nodes by Y
 *      K : number of nodes by time
 *      x_0, y_0 : X and Y values at left boundaries
 *      x_1, y_1 : X and Y values at right boundaries
 *      t_1 : end time value
 *      hx, hy, tau : step by X, Y and time respectively
 */
extern int N;
extern int M;
extern int K;
float x_scale = 1;
float y_scale = 1;
float time_scale = 1;
float x_0 = 0.0;
float y_0 = 0.0;
float x_1 = 1.0;
float y_1 = 1.0;
float t_1 = 1.0;
float hx;
float hy;
float tau;
float T0 = 1.0/4.0;
float T1 = 1.0;


/*
 * Defining constants:
 *      Re : Reynolds number
 *      Pr : Prandtl number
 *      Gr : Grashof number
 *      xi : coefficient of thermal conductivity
 *      nu : coefficient of kinematic viscosity
 *      right_beta : coefficient of thermal expansion
 */
float Re = 1.0;
float Gr = 10000.0;
float Pr = 1.0;
float xi;
float nu;
float rigth_beta;

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

void init_boundary_conditions(Matrix3D array, char side, int a, int b, float value) {
    switch (side)
    {
    case 'l':
        for (int k = 0; k < K; k++)
            for (int j = 0; j < M; j++)
                array[k][j][0] = a * array[k][j][1] + b*value;
        break;

    case 'r':
        for (int k = 0; k < K; k++)
            for (int j = 0; j < M; j++)
                array[k][j][N-1] = a * array[k][j][N-2] + b*value;
        break;
    
    case 'b':
        for (int k = 0; k < K; k++)
            for (int i = 0; i < N; i++)
                array[k][0][i] = a * array[k][1][i] + b*value;
        break;
    
    case 't':
        for (int k = 0; k < K; k++)
            for (int i = 0; i < N; i++)
                array[k][M-1][i] = a * array[k][M-2][i] + b*value;
        break;

    default:
        break;
    }
}

void calculations() {
    /*
     * Global cycle for time
     */
    for (int k = 0; k < K-1; k++)
    {
        /*
         * Run back-for by X for over Y on n+1/2 layer
         */
        float coef_x = tau/(2*hx);
        for (int j = 1; j < M-1; j++)
        {
            alpha_x[1] = 0;
            beta_x[1] = T0;

            for (int i = 1; i < N-1; i++)
            {
                float A = coef_x * (xi/hx + u[j][i]/2);
                float B = coef_x * (xi/hx - u[j][i]/2);
                float C = 1 + tau * xi / (hx * hx);
                float F = temperature[k][j][i] - tau*v[j][i]/(4*hy)*(temperature[k][j+1][i] - temperature[k][j-1][i]) + \
                            tau*xi/(2*hy*hy)*(temperature[k][j+1][i] - 2*temperature[k][j][i] + temperature[k][j-1][i]);
                

                alpha_x[i+1] = B / (C - A * alpha_x[i]);
                beta_x[i+1] = (A * beta_x[i] + F) / (C - A * alpha_x[i]);
            }

            temperature12[k][j][N-1] = T1;
            for (int i = N-2; i >= 0; i--)
            {
                temperature12[k][j][i] = alpha_x[i+1] * temperature12[k][j][i+1] + beta_x[i+1];
            }
        }

        /*
         * Run back-for by Y for over X on n+1 layer
         */
        float coef_y = tau/(2*hy);
        for (int i = 1; i < N-1; i++)
        {
            // alpha - a second-order condition.
            // beta - a first-order condition.
            alpha_y[1] = 1;
            beta_y[1] = 0;

            for (int j = 1; j < M-1; j++)
            {
                float A = coef_y * (xi/hy + v[j][i]/2);
                float B = coef_y * (xi/hy - v[j][i]/2);
                float C = 1 + tau * xi / (hy * hy);
                float F = temperature12[k][j][i] - tau*u[j][i]/(4*hx)*(temperature12[k][j][i+1] - temperature12[k][j][i-1]) + \
                            tau*xi/(2*hx*hx)*(temperature12[k][j][i+1] - 2*temperature12[k][j][i] + temperature12[k][j][i-1]);

                alpha_y[j+1] = B / (C - A * alpha_y[j]);
                beta_y[j+1] = (A * beta_y[j] + F) / (C - A * alpha_y[j]);
            }

            temperature[k+1][M-1][i] = beta_y[N-1]/(1-alpha_y[N-1]);
            for (int j = M-2; j >= 0; j--)
            {
                temperature[k+1][j][i] = alpha_y[j+1] * temperature[k+1][j+1][i] + beta_y[j+1];
                // temperature[k+1][j][0] = temperature[k+1][j][1];
                // temperature[k+1][j][N-1] = temperature[k+1][j][N-2];
            }
        }
        
        /*
         * Run back-for by X for over Y on n+1/2 layer
         */
        for (int j = 1; j < M-1; j++)
        {
            // alpha - a second-order condition.
            // beta - a first-order condition.
            alpha_x[1] = 0;
            beta_x[1] = -2 / (hx*hx) * psi[j][1];

            for (int i = 1; i < N-1; i++)
            {
                float A = coef_x * (nu/hx + u[j][i]/2);
                float B = coef_x * (nu/hx - u[j][i]/2);
                float C = 1 + tau * nu / (hx * hx);
                float F = omega[k][j][i] - tau*v[j][i]/(4*hy)*(omega[k][j+1][i] - omega[k][j-1][i]) + \
                            tau*nu/(2*hy*hy)*(omega[k][j+1][i] - 2*omega[k][j][i] + omega[k][j-1][i]) + \
                            tau*rigth_beta/(4*hx*hx)*(temperature[k][j][i+1] - temperature[k][j][i-1]);

                alpha_x[i+1] = B / (C - A * alpha_x[i]);
                beta_x[i+1] = (A * beta_x[i] + F) / (C - A * alpha_x[i]);
            }

            omega12[k][j][N-1] = -2 / (hx*hx) * psi[j][N-2];
            for (int i = N-2; i >= 0; i--)
                omega12[k][j][i] = alpha_x[i+1] * omega12[k][j][i+1] + beta_x[i+1];
        }

        /*
         * Run back-for by Y for over X on n+1 layer
         */
        for (int i = 1; i < N-1; i++)
        {
            // alpha - a second-order condition.
            // beta - a first-order condition.
            alpha_y[1] = 0;
            beta_y[1] = -2 / (hy*hy) * psi[1][i];

            for (int j = 1; j < M-1; j++)
            {
                float A = coef_y * (nu/hy - v[j][i]/2);
                float B = coef_y * (nu/hy + v[j][i]/2);
                float C = 1 + tau * nu / (hy * hy);
                float F = omega12[k][j][i] - tau*u[j][i]/(4*hx)*(omega12[k][j][i+1] - omega12[k][j][i-1]) + \
                            tau*nu/(2*hx*hx)*(omega12[k][j][i+1] - 2*omega12[k][j][i] + omega12[k][j][i-1]) + \
                            tau*rigth_beta/(4*hx*hx)*(temperature[k][j][i+1] - temperature[k][j][i-1]);

                alpha_y[j+1] = B / (C - A * alpha_y[j]);
                beta_y[j+1] = (A * beta_y[j] + F) / (C - A * alpha_y[j]);
            }

            omega[k+1][M-1][i] = -2 / (hy*hy) * psi[M-2][i];
            for (int j = M-2; j >= 0; j--)
            {
                omega[k+1][j][i] = alpha_y[j+1] * omega[k+1][j+1][i] + beta_y[j+1];
                omega[k+1][j][N-1] = -2 / (hy*hy) * psi[j][N-2];
                omega[k+1][j][0] = -2 / (hy*hy) * psi[j][1];
            }
        }

        /*
         * Calculations of PSI with star-scheme
         */
        float psi_coef = (hx*hx*hy*hy/(2*hx*hx + 2*hy*hy));
        for (int j = 1; j < M-1; j++)
        {
            for (int i = 1; i < N-1; i++)
            {
                psi[j][i] = psi_coef * (
                                (psi[j][i+1] + psi[j][i-1])/(hx*hx) + \
                                (psi[j+1][i] + psi[j-1][i])/(hy*hy) + \
                                omega[k][j][i]);
            }
        }

        /*
         * Calculations of velocities
         */
        for (int j = 1; j < M-1; j++)
        {
            for (int i = 1; i < N-1; i++)
            {
                u[j][i] = (psi[j+1][i] - psi[j-1][i]) / (2 * hy);
                v[j][i] = -(psi[j][i+1] - psi[j][i-1]) / (2 * hx);
            }
        }
    }
}

void helper(const char *progname) {
    printf("%s starts a calculations convective flow of fluid in 2D case.\n\n", progname);
    printf("Usage:\n");
    printf("    %s [OPTION]...      \n", progname);
    printf("\nOptions:\n");
    printf("    -G, --grashof       Reynolds number\n");
    printf("    -P, --prandtl       Prandtl number\n");
    printf("    -R, --reynolds      Grashof number\n");
    printf("    -T, --time-points   number of nodes by X\n");
    printf("    -t, --t1            \n");
    printf("    -x, --x1            \n");
    printf("    -X, --x-points      number of nodes by Y\n");
    printf("    -y, --y1            \n");
    printf("    -Y, --y-points      number of nodes by time\n");
    printf("        --x-scale       \n");
    printf("        --y-scale       \n");
    printf("    -?, --help          \n");
}

void free_arrays(){
    printf("\nFREEE ALL ALLOCATED MEMORY\n");
    free_vector(x);
    free_vector(y);
    free_vector(t);
    free_matrix3d(temperature);
    free_matrix3d(temperature12);
    free_matrix3d(omega);
    free_matrix3d(omega12);
    free_matrix(psi);
    free_matrix(u);
    free_matrix(v);
    free_vector(alpha_x);
    free_vector(beta_x);
    free_vector(alpha_y);
    free_vector(beta_y);
    
}

int main(int argc, char *argv[]) {

    if (argc < 1)
    {
        printf("No arguments are given!");
        return 0;
    }

    static struct option long_options[] = {
        {"grashof", required_argument, NULL, 'G'},
        {"help", no_argument, NULL, '?'},
        {"prandtl", required_argument, NULL, 'P'},
        {"reynolds", required_argument, NULL, 'R'},
        {"t_1", required_argument, NULL, 't'},
        {"time-points", required_argument, NULL, 'T'},
        {"x_1", required_argument, NULL, 'x'},
        {"x-points", required_argument, NULL, 'X'},
        {"x-scale", required_argument, NULL, 1},
        {"y_1", required_argument, NULL, 'y'},
        {"y-points", required_argument, NULL, 'Y'},
        {"y-scale", required_argument, NULL, 2},
        {NULL, 0, NULL, 0}
	};
    int c;
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "G:P:R:t:T:x:X:y:Y:",
                 long_options, &option_index)) != -1)
    {
        switch (c)
        {
            case 1:
                x_scale = atof(optarg);
                NUMCHECK(x_scale)
                break;
            case 2:
                y_scale = atof(optarg);
                NUMCHECK(y_scale)
                break;
            case 'X':
                N = atoi(optarg);
                if (N < 10)
                {
                    printf("Wrong number of X steps: %d\n", N);
                    exit(1);
                }
                break;
            case 'Y':
                M = atoi(optarg);
                if (M < 10)
                {
                    printf("Wrong number of Y steps: %d\n", M);
                    exit(1);
                }
                break;
            case 'T':
                K = atoi(optarg);
                if (K < 10)
                {
                    printf("Wrong number of time steps: %d\n", K);
                    exit(1);
                }
                break;
            case 'x':
                x_1 = atof(optarg);
                NUMCHECK(x_1)
                break;
            case 'y':
                y_1 = atof(optarg);
                NUMCHECK(y_1)
                break;
            case 't':
                t_1 = atof(optarg);
                NUMCHECK(t_1)
                break;
            case 'R':
                Re = atof(optarg);
                NUMCHECK(Re)
                break;
            case 'G':
                Gr = atof(optarg);
                NUMCHECK(Gr)
                break;
            case 'P':
                Pr = atof(optarg);
                NUMCHECK(Pr)
                break;
            case '?':
                helper(argv[0]);
                exit(0);
                break;
            default:
                printf("Try \"./convective_flow --help\" for more information.\n");
                exit(1);
        }
    }

    hx = x_1 / N;
    hy = y_1 / M;
    tau = t_1 / K;
    xi = 1.0 / (Re * Pr);
    nu = 1.0 / Re;
    rigth_beta = Gr / (Re * Re);

    x = init_zero_vector(N+1);
    y = init_zero_vector(M+1);
    t = init_zero_vector(K);
    temperature = init_zero_matrix3d(N, M, K);
    temperature12 = init_zero_matrix3d(N, M, K);
    omega = init_zero_matrix3d(N, M, K);
    omega12 = init_zero_matrix3d(N, M, K);
    psi = init_zero_matrix(N, M);
    u = init_zero_matrix(N, M);
    v = init_zero_matrix(N, M);
    alpha_x = init_zero_vector(N);
    beta_x = init_zero_vector(N);
    alpha_y = init_zero_vector(M);
    beta_y = init_zero_vector(M);
    if ( x == NULL || y == NULL || t == NULL || \
         temperature == NULL || temperature12 == NULL || \
        omega == NULL || omega12 == NULL || psi == NULL || \
        u == NULL || v == NULL || alpha_x == NULL || beta_y == NULL )
    {
        printf("Failed to allocate memory\n");
        free_arrays();
        exit(1);
    }

    for (int i = 0; i < N+1; i++)
        x[i] = i*hx;

    for (int j = 0; j < M+1; j++)
        y[j] = j*hy;

    for (int k = 0; k < K; k++)
        t[k] = k*tau;

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            temperature[0][j][i] = exp(-190*(x[i])) + T0;
            temperature12[0][j][i] = exp(-190*(x[i])) + T0;
        }
    }

    init_boundary_conditions(temperature, 'l', 0, 1, T0);
    init_boundary_conditions(temperature, 'r', 0, 1, T1);
    init_boundary_conditions(temperature12, 'l', 0, 1, T0);
    init_boundary_conditions(temperature12, 'r', 0, 1, T1);

    printf("RUN CALCULATIONS\n");
    calculations();
    printf("\tSuccessful!\n");
    
    printf("SAVE DATA\n");
    FILE* temp = fopen("./output_data/Temp.txt", "w");
    FILE* Omg = fopen("./output_data/Omg.txt", "w");
    FILE* Psi = fopen("./output_data/Psi.txt", "w");
    if (temp == NULL || Omg == NULL || Psi == NULL)
        return -1;

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            fprintf(temp, "%.4f ", temperature[K-1][j][i]);
            fprintf(Omg, "%.4f ", omega[K-1][j][i]);
            fprintf(Psi, "%.4f ", psi[j][i]);
        }
        fprintf(temp, "%c", '\n');
        fprintf(Omg, "%c", '\n');
        fprintf(Psi, "%c", '\n');
    }

    fclose(temp);
    fclose(Omg);
    fclose(Psi);
    printf("\tSuccessful!\n");

    printf("SAVE TIME DATA\n");
    FILE* temp_all = fopen("./output_data/temp_all.txt", "w");
    if (temp_all == NULL)
        return -1;
    
    for (int k = 0; k < K; k+=10)
    {
        for (int j = 0; j < M; j++)
        {
            for (int i = 0; i < N; i++)
            {
                fprintf(temp_all, "%.4f ", temperature[k][j][i]);
            }
        }
        fprintf(temp_all, "%c", '\n');
    }

    fclose(temp_all);
    printf("\tSuccessful!\n");
    free_arrays();
    printf("\tSuccessful!\n");

    return 0;
}
