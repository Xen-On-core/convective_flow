#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include "utils/array.h"
#include "utils/grid.h"
#include "utils/utils.h"
#include "utils/solutions.h"



double epsilon = 1e-6;
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
int N;
int M;
int K;
int T;
double x_0 = 0.0;
double y_0 = 0.0;
double x_1 = 1.0;
double y_1 = 1.0;
double t_1 = 1.0;
double hx;
double hy;
double tau;
double T0 = 1.0/4.0;
double T1 = 1.0;

CFP cfparams;
Grid2D grid;
CFS_2D solution;

Vector t;
Matrix3D temperature;
Matrix3D temperature12;
Matrix3D omega;
Matrix3D omega12;
Matrix3D psi;
Matrix3D psi12;
Matrix u;
Matrix v;
Vector alpha_x;
Vector beta_x;
Vector alpha_y;
Vector beta_y;

void init_boundary_conditions(Matrix3D array, char side, int a, int b, double value)
{
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

void calculations()
{
    /*
     * Global cycle for time
     */
    for (int k = 0; k < K-1; k++)
    {
        /*
         * Run back-for by X for over Y on n+1/2 layer
         */
        double coef_x = tau/(2*hx);
        for (int j = 1; j < M-1; j++)
        {
            alpha_x[1] = 0;
            beta_x[1] = T0;

            for (int i = 1; i < N-1; i++)
            {
                double A = coef_x * (cfparams.xi/hx + u[j][i]/2);
                double B = coef_x * (cfparams.xi/hx - u[j][i]/2);
                double C = 1 + tau * cfparams.xi / (hx * hx);
                double F = temperature[k][j][i] - tau*v[j][i]/(4*hy)*(temperature[k][j+1][i] - temperature[k][j-1][i]) + \
                            tau*cfparams.xi/(2*hy*hy)*(temperature[k][j+1][i] - 2*temperature[k][j][i] + temperature[k][j-1][i]);
                

                alpha_x[i+1] = B / (C - A * alpha_x[i]);
                beta_x[i+1] = (A * beta_x[i] + F) / (C - A * alpha_x[i]);
            }

            temperature12[k][j][N-1] = T1; //(tau*k >= 0.5) ? (T1+1.0) : T1 + 2*tau*k;
            for (int i = N-2; i >= 0; i--)
            {
                temperature12[k][j][i] = alpha_x[i+1] * temperature12[k][j][i+1] + beta_x[i+1];
            }
        }

        /*
         * Run back-for by Y for over X on n+1 layer
         */
        double coef_y = tau/(2*hy);
        for (int i = 1; i < N-1; i++)
        {
            // alpha - a second-order condition.
            // beta - a first-order condition.
            alpha_y[1] = 1;
            beta_y[1] = 0;

            for (int j = 1; j < M-1; j++)
            {
                double A = coef_y * (cfparams.xi/hy + v[j][i]/2);
                double B = coef_y * (cfparams.xi/hy - v[j][i]/2);
                double C = 1 + tau * cfparams.xi / (hy * hy);
                double F = temperature12[k][j][i] - tau*u[j][i]/(4*hx)*(temperature12[k][j][i+1] - temperature12[k][j][i-1]) + \
                            tau*cfparams.xi/(2*hx*hx)*(temperature12[k][j][i+1] - 2*temperature12[k][j][i] + temperature12[k][j][i-1]);

                alpha_y[j+1] = B / (C - A * alpha_y[j]);
                beta_y[j+1] = (A * beta_y[j] + F) / (C - A * alpha_y[j]);
            }

            temperature[k+1][M-1][i] = beta_y[N-1]/(1-alpha_y[N-1]);
            for (int j = M-2; j >= 0; j--)
            {
                temperature[k+1][j][i] = alpha_y[j+1] * temperature[k+1][j+1][i] + beta_y[j+1];
                // temperature[k+1][j][0] = T0;
                // temperature[k+1][j][N-1] = (k >= 0.5) ? (T1+1) : (T1 + 2*tau*(k+1));
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
            beta_x[1] = -2 / (hx*hx) * psi[k][j][1];

            for (int i = 1; i < N-1; i++)
            {
                double A = coef_x * (cfparams.nu/hx + u[j][i]/2);
                double B = coef_x * (cfparams.nu/hx - u[j][i]/2);
                double C = 1 + tau * cfparams.nu / (hx * hx);
                double F = omega[k][j][i] - tau*v[j][i]/(4*hy)*(omega[k][j+1][i] - omega[k][j-1][i]) + \
                            tau*cfparams.nu/(2*hy*hy)*(omega[k][j+1][i] - 2*omega[k][j][i] + omega[k][j-1][i]) + \
                            tau*cfparams.rigth_beta/(4*hx*hx)*(temperature[k][j][i+1] - temperature[k][j][i-1]);

                alpha_x[i+1] = B / (C - A * alpha_x[i]);
                beta_x[i+1] = (A * beta_x[i] + F) / (C - A * alpha_x[i]);
            }

            omega12[k][j][N-1] = -2 / (hx*hx) * psi[k][j][N-2];
            for (int i = N-2; i >= 0; i--)
            {
                omega12[k][j][i] = alpha_x[i+1] * omega12[k][j][i+1] + beta_x[i+1];
            }
            
        }
        /*
         * Run back-for by Y for over X on n+1 layer
         */
        for (int i = 1; i < N-1; i++)
        {
            // alpha - a second-order condition.
            // beta - a first-order condition.
            alpha_y[1] = 0;
            beta_y[1] = -2 / (hy*hy) * psi[k][1][i];

            for (int j = 1; j < M-1; j++)
            {
                double A = coef_y * (cfparams.nu/hy - v[j][i]/2);
                double B = coef_y * (cfparams.nu/hy + v[j][i]/2);
                double C = 1 + tau * cfparams.nu / (hy * hy);
                double F = omega12[k][j][i] - tau*u[j][i]/(4*hx)*(omega12[k][j][i+1] - omega12[k][j][i-1]) + \
                            tau*cfparams.nu/(2*hx*hx)*(omega12[k][j][i+1] - 2*omega12[k][j][i] + omega12[k][j][i-1]) + \
                            tau*cfparams.rigth_beta/(4*hx*hx)*(temperature12[k][j][i+1] - temperature12[k][j][i-1]);

                alpha_y[j+1] = B / (C - A * alpha_y[j]);
                beta_y[j+1] = (A * beta_y[j] + F) / (C - A * alpha_y[j]);
            }

            omega[k+1][M-1][i] = -2 / (hy*hy) * psi[k][M-2][i];
            for (int j = M-2; j >= 0; j--)
            {
                omega[k+1][j][i] = alpha_y[j+1] * omega[k+1][j+1][i] + beta_y[j+1];
                omega[k+1][j][N-1] = -2 / (hy*hy) * psi[k][j][N-2];
                omega[k+1][j][0] = -2 / (hy*hy) * psi[k][j][1];
            }
        }

        /*
         * Calculations of PSI
         */
        for (int j = 1; j < M-1; j++)
        {
	        alpha_x[1] = 1;
	        beta_x[1] = 0;
            
            for (int i = 1; i < N-1; i++)
            {
                double A = tau * cfparams.lambda/(2*hx*hx);
                double B = tau * cfparams.lambda/(2*hx*hx);
		        double C = 1 + A + B;
		        double F = psi[k][j][i] + tau*cfparams.lambda/(2*hy*hy)*(psi[k][j+1][i] - 2*psi[k][j][i] + psi[k][j-1][i]) + tau*cfparams.lambda/2 * omega[k][j][i];

		        alpha_x[i+1] = B / (C - alpha_x[i] * A);
		        beta_x[i+1] = (F + A * beta_x[i]) / (C - alpha_x[i] * A);
            }

            psi12[k][j][N-1] = beta_x[N-1]/(1-alpha_x[N-1]);
            for (int i = N-2; i >= 0; i--)
            {
                psi12[k][j][i] = alpha_x[i+1] * psi12[k][j][i+1] + beta_x[i+1];
            }
        }

	    for (int i=1; i < N-1; i++)
	    {
	        alpha_y[1] = 1;
	        beta_y[1] = 0;

            for (int j = 1; j < M-1; j++)
            {
                double A = tau * cfparams.lambda/(2*hy*hy);
                double B = tau * cfparams.lambda/(2*hy*hy);
	    	    double C = 1 + A + B;
	    	    double F = psi12[k][j][i] + tau*cfparams.lambda/(2*hx*hx)*(psi12[k][j][i+1] - 2*psi12[k][j][i] + psi12[k][j][i-1]) + tau*cfparams.lambda/2 * omega[k][j][i];

	    	    alpha_y[j+1] = B / (C - alpha_y[j] * A);
	    	    beta_y[j+1] = (F + A * beta_y[j]) / (C - alpha_y[j] * A);
            }

            psi[k+1][M-1][i] = beta_y[M-1]/(1-alpha_y[M-1]);
            for (int j = M-2; j >= 0; j--)
            {
                psi[k+1][j][i] = alpha_y[j+1] * psi[k+1][j+1][i] + beta_y[j+1];
            }
	    }


        /*
         * Calculations of velocities
         */
        for (int j = 1; j < M-1; j++)
        {
            for (int i = 1; i < N-1; i++)
            {
                u[j][i] = (psi[k][j+1][i] - psi[k][j-1][i]) / (2 * hy);
                v[j][i] = -(psi[k][j][i+1] - psi[k][j][i-1]) / (2 * hx);
            }
        }
    }
}

void helper(const char *progname)
{
    printf("%s starts a calculations convective flow of fluid in 2D case.\n\n", progname);
    printf("Usage:\n");
    printf("    %s [OPTION]...      \n", progname);
    printf("\nOptions:\n");
    printf("    -G, --grashof       Reynolds number\n");
    printf("    -P, --prandtl       Prandtl number\n");
    printf("    -R, --reynolds      Grashof number\n");
    printf("    -T, --time-points   number of nodes by time\n");
    printf("    -t, --t1            \n");
    printf("    -x, --x1            \n");
    printf("    -X, --x-points      number of nodes by X\n");
    printf("    -y, --y1            \n");
    printf("    -Y, --y-points      number of nodes by Y\n");
    printf("        --x-scale       \n");
    printf("        --y-scale       \n");
    printf("    -?, --help          \n");
}

void free_arrays()
{
    printf("\nFREEE ALL ALLOCATED MEMORY\n");
    free_vector(grid.x1);
    free_vector(grid.x2);
    free_vector(t);
    free_matrix3d(temperature);
    free_matrix3d(temperature12);
    free_matrix3d(omega);
    free_matrix3d(omega12);
    free_matrix3d(psi);
    free_matrix3d(psi12);
    free_matrix(u);
    free_matrix(v);
    free_vector(alpha_x);
    free_vector(beta_x);
    free_vector(alpha_y);
    free_vector(beta_y);
}

void save_data()
{
    printf("SAVE DATA\n");
    FILE* temperature_file = fopen("./output_data/temperature.txt", "w");
    FILE* omega_file = fopen("./output_data/omega.txt", "w");
    FILE* psi_file = fopen("./output_data/psi.txt", "w");
    FILE* velocity_u = fopen("./output_data/u_vel.txt", "w");
    FILE* velocity_v = fopen("./output_data/v_vel.txt", "w");
    if (temperature == NULL || omega_file == NULL || psi_file == NULL)
        return;

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            fprintf(temperature_file, "%.4f ", temperature[K-1][j][i]);
            fprintf(omega_file, "%.4f ", omega[K-1][j][i]);
            fprintf(psi_file, "%.4f ", psi[K-1][j][i]);
            fprintf(velocity_u, "%.4f ", u[j][i]);
            fprintf(velocity_v, "%.4f ", v[j][i]);
        }
        fprintf(temperature_file, "%c", '\n');
        fprintf(omega_file, "%c", '\n');
        fprintf(psi_file, "%c", '\n');
        fprintf(velocity_u, "%c", '\n');
        fprintf(velocity_v, "%c", '\n');
    }

    fclose(temperature_file);
    fclose(omega_file);
    fclose(psi_file);
    fclose(velocity_u);
    fclose(velocity_v);
    printf("\tSuccessful!\n");

    printf("SAVE TIME DATA\n");
    FILE* temperature_time_data = fopen("./output_data/temperature_time_data.txt", "w");
    if (temperature_time_data == NULL)
        return;
    
    for (int k = 0; k < K; k+=10)
    {
        for (int j = 0; j < M; j++)
        {
            for (int i = 0; i < N; i++)
            {
                fprintf(temperature_time_data, "%.4f ", temperature[k][j][i]);
            }
        }
        fprintf(temperature_time_data, "%c", '\n');
    }
    fclose(temperature_time_data);
}

double getMaxAbs(Matrix3D p)
{
    Matrix subtr = init_zero_matrix(N, M);

    double max = 0.0;
    double maxitem = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            subtr[j][i] = fabs(p[K-1][j][i] - p[K-2][j][i]);
            if (max <= subtr[j][i])
                max = subtr[j][i];

	    if (maxitem <= fabs(p[K-1][j][i]))
	        maxitem = p[K-1][j][i];
        }
    }

    printf("MAX PSI ITEM: %f\n", maxitem);
    
    return max;
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
        {"y_1", required_argument, NULL, 'y'},
        {"y-points", required_argument, NULL, 'Y'},
        {NULL, 0, NULL, 0}
	};
    int c;
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "G:P:R:t:T:x:X:y:Y:",
                 long_options, &option_index)) != -1)
    {
        switch (c)
        {
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
                grid.x1_r = atof(optarg);
                NUMCHECK(grid.x1_r);
                break;
            case 'y':
                grid.x2_r = atof(optarg);
                NUMCHECK(grid.x2_r);
                break;
            case 't':
                t_1 = atof(optarg);
                NUMCHECK(t_1);
                break;
            case 'R':
                cfparams.Re = atof(optarg);
                NUMCHECK(cfparams.Re);
                break;
            case 'G':
                cfparams.Gr = atof(optarg);
                NUMCHECK(cfparams.Gr);
                break;
            case 'P':
                cfparams.Pr = atof(optarg);
                NUMCHECK(cfparams.Pr);
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

    K = K * t_1 + 1;
    grid.x1_size = N * x_1 + 1;
    grid.x2_size = M * y_1 + 1;

    tau = t_1 / (K-1);
    grid.hx1 = x_1 / (grid.x1_size-1);
    grid.hx2 = y_1 / (grid.x2_size-1);
    hx = grid.hx1;
    hy = grid.hx2;

    cfparams.xi = 1.0 / (cfparams.Re * cfparams.Pr);
    cfparams.nu = 1.0 / cfparams.Re;
    cfparams.rigth_beta = cfparams.Gr / (cfparams.Re * cfparams.Re);
    cfparams.lambda = 0.5;

    t = init_zero_vector(K);
    grid.x1 = init_zero_vector(grid.x1_size);
    grid.x2 = init_zero_vector(grid.x2_size);
    alpha_x = init_zero_vector(grid.x1_size);
    alpha_y = init_zero_vector(grid.x2_size);
    beta_x = init_zero_vector(grid.x1_size);
    beta_y = init_zero_vector(grid.x2_size);
    u = init_zero_matrix(grid.x1_size, grid.x2_size);
    v = init_zero_matrix(grid.x1_size, grid.x2_size);
    psi = init_zero_matrix3d(grid.x1_size, grid.x2_size, K);
    psi12 = init_zero_matrix3d(grid.x1_size, grid.x2_size, K);
    omega = init_zero_matrix3d(grid.x1_size, grid.x2_size, K);
    omega12 = init_zero_matrix3d(grid.x1_size, grid.x2_size, K);
    temperature = init_zero_matrix3d(grid.x1_size, grid.x2_size, K);
    temperature12 = init_zero_matrix3d(grid.x1_size, grid.x2_size, K);
    if ( grid.x1 == NULL || grid.x2 == NULL || t == NULL || \
         temperature == NULL || temperature12 == NULL || \
        omega == NULL || omega12 == NULL || psi == NULL || \
        u == NULL || v == NULL || alpha_x == NULL || beta_y == NULL )
    {
        printf("Failed to allocate memory\n");
        free_arrays();
        exit(1);
    }

    for (int i = 0; i < grid.x1_size; i++)
        grid.x1[i] = i*hx;

    for (int j = 0; j < grid.x2_size; j++)
        grid.x2[j] = j*hy;

    for (int k = 0; k < K; k++)
        t[k] = k*tau;

    // for (int j = 0; j < M; j++)
    // {
    //     for (int i = 0; i < N; i++)
    //     {
    //         temperature[0][j][i] = exp(-190*(x[i])) + T0;
    //         temperature12[0][j][i] = exp(-190*(x[i])) + T0;
    //     }
    // }

    init_boundary_conditions(temperature, 'l', 0, 1, T0);
    init_boundary_conditions(temperature, 'r', 0, 1, T1);
    init_boundary_conditions(temperature12, 'l', 0, 1, T0);
    init_boundary_conditions(temperature12, 'r', 0, 1, T1);

    printf("RUN CALCULATIONS\n");
    calculations();
    print_matrix(psi[5]);
    printf("\tSuccessful!\n");
    
    save_data();
    
    printf("%f\n", getMaxAbs(psi));
    printf("\tSuccessful!\n");
    free_arrays();
    printf("\tSuccessful!\n");

    return 0;
}
