#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include "utils/array.h"
#include "utils/grid.h"
#include "utils/utils.h"

FILE* logfile;

double epsilon = 1e-6;
#define NUMCHECK(number)                                              \
    if (abs(number) < epsilon)                                        \
    {                                                                 \
        fprintf(logfile, "Wrong value of %s: %f\n", #number, number); \
        exit(1);                                                      \
    }

int alloc_error = 0;

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
int T;
int K;
int N;
int M;
double tau;
double T0 = 1.0 / 4.0;
double T1 = 1.0;

CFP cfparams;
Grid2D grid;

Vector *t;
Matrix3D *temperature;
Matrix3D *temperature12;
Matrix3D *omega;
Matrix3D *omega12;
Matrix3D *psi;
Matrix3D *psi12;
Matrix *u;
Matrix *v;
Vector *alpha_x;
Vector *beta_x;
Vector *alpha_y;
Vector *beta_y;

void calculations()
{
    for (int k = 0; k < K - 1; k++)
    {
        double coef_x = tau / (2 * grid.hx1);
        for (int j = 1; j < grid.x1->size - 1; j++)
        {
            alpha_x->data[1] = 0;
            beta_x->data[1] = T0;

            for (int i = 1; i < grid.x1->size - 1; i++)
            {
                double A = coef_x * (cfparams.xi / grid.hx1 + u->data[j][i] / 2);
                double B = coef_x * (cfparams.xi / grid.hx1 - u->data[j][i] / 2);
                double C = 1 + tau * cfparams.xi / (grid.hx1 * grid.hx1);
                double F = temperature->data[k][j][i] - tau * v->data[j][i] / (4 * grid.hx2) * 
                            (temperature->data[k][j + 1][i] - temperature->data[k][j - 1][i]) + 
                            tau * cfparams.xi / (2 * grid.hx2 * grid.hx2) * 
                            (temperature->data[k][j + 1][i] - 2 * temperature->data[k][j][i] + temperature->data[k][j - 1][i]);

                alpha_x->data[i + 1] = B / (C - A * alpha_x->data[i]);
                beta_x->data[i + 1] = (A * beta_x->data[i] + F) / (C - A * alpha_x->data[i]);
            }

            temperature12->data[k][j][grid.x1->size - 1] = T1; //(tau*k >= 0.5) ? (T1+1.0) : T1 + 2*tau*k;
            for (int i = grid.x1->size - 2; i >= 0; i--)
            {
                temperature12->data[k][j][i] = alpha_x->data[i + 1] * temperature12->data[k][j][i + 1] + beta_x->data[i + 1];
            }
        }

        double coef_y = tau / (2 * grid.hx2);
        for (int i = 1; i < grid.x1->size - 1; i++)
        {
            alpha_y->data[1] = 1;
            beta_y->data[1] = 0;

            for (int j = 1; j < grid.x2->size - 1; j++)
            {
                double A = coef_y * (cfparams.xi / grid.hx2 + v->data[j][i] / 2);
                double B = coef_y * (cfparams.xi / grid.hx2 - v->data[j][i] / 2);
                double C = 1 + tau * cfparams.xi / (grid.hx2 * grid.hx2);
                double F = temperature12->data[k][j][i] - tau * u->data[j][i] / (4 * grid.hx1) * 
                            (temperature12->data[k][j][i + 1] - temperature12->data[k][j][i - 1]) + 
                            tau * cfparams.xi / (2 * grid.hx1 * grid.hx1) * 
                            (temperature12->data[k][j][i + 1] - 2 * temperature12->data[k][j][i] + temperature12->data[k][j][i - 1]);

                alpha_y->data[j + 1] = B / (C - A * alpha_y->data[j]);
                beta_y->data[j + 1] = (A * beta_y->data[j] + F) / (C - A * alpha_y->data[j]);
            }

            temperature->data[k + 1][grid.x2->size - 1][i] = beta_y->data[grid.x2->size - 1] / (1 - alpha_y->data[grid.x2->size - 1]);
            temperature->data[k + 1][grid.x2->size - 1][0] = T0;
            temperature->data[k + 1][grid.x2->size - 1][grid.x1->size - 1] = T1;
            for (int j = grid.x2->size - 2; j >= 0; j--)
            {
                temperature->data[k + 1][j][i] = alpha_y->data[j + 1] * temperature->data[k + 1][j + 1][i] + beta_y->data[j + 1];
                temperature->data[k + 1][j][0] = T0;
                temperature->data[k + 1][j][grid.x1->size - 1] = T1; //(k >= 0.5) ? (T1+1) : (T1 + 2*tau*(k+1));
            }
        }

        for (int j = 1; j < grid.x2->size - 1; j++)
        {
            alpha_x->data[1] = 0;
            beta_x->data[1] = -2 / (grid.hx1 * grid.hx1) * psi->data[k][j][1];

            for (int i = 1; i < grid.x1->size - 1; i++)
            {
                double A = coef_x * (cfparams.nu / grid.hx1 + u->data[j][i] / 2);
                double B = coef_x * (cfparams.nu / grid.hx1 - u->data[j][i] / 2);
                double C = 1 + tau * cfparams.nu / (grid.hx1 * grid.hx1);
                double F = omega->data[k][j][i] - tau * v->data[j][i] / (4 * grid.hx2) * 
                            (omega->data[k][j + 1][i] - omega->data[k][j - 1][i]) + 
                            tau * cfparams.nu / (2 * grid.hx2 * grid.hx2) * 
                            (omega->data[k][j + 1][i] - 2 * omega->data[k][j][i] + omega->data[k][j - 1][i]) + 
                            tau * cfparams.rigth_beta / (4 * grid.hx1 * grid.hx1) * 
                            (temperature->data[k][j][i + 1] - temperature->data[k][j][i - 1]);

                alpha_x->data[i + 1] = B / (C - A * alpha_x->data[i]);
                beta_x->data[i + 1] = (A * beta_x->data[i] + F) / (C - A * alpha_x->data[i]);
            }

            omega12->data[k][j][grid.x1->size - 1] = -2 / (grid.hx1 * grid.hx1) * psi->data[k][j][grid.x1->size - 2];
            for (int i = grid.x1->size - 2; i >= 0; i--)
            {
                omega12->data[k][j][i] = alpha_x->data[i + 1] * omega12->data[k][j][i + 1] + beta_x->data[i + 1];
            }
        }

        for (int i = 1; i < grid.x1->size - 1; i++)
        {
            alpha_y->data[1] = 0;
            beta_y->data[1] = -2 / (grid.hx2 * grid.hx2) * psi->data[k][1][i];

            for (int j = 1; j < grid.x2->size - 1; j++)
            {
                double A = coef_y * (cfparams.nu / grid.hx2 - v->data[j][i] / 2);
                double B = coef_y * (cfparams.nu / grid.hx2 + v->data[j][i] / 2);
                double C = 1 + tau * cfparams.nu / (grid.hx2 * grid.hx2);
                double F = omega12->data[k][j][i] - tau * u->data[j][i] / (4 * grid.hx1) * 
                            (omega12->data[k][j][i + 1] - omega12->data[k][j][i - 1]) + 
                            tau * cfparams.nu / (2 * grid.hx1 * grid.hx1) * 
                            (omega12->data[k][j][i + 1] - 2 * omega12->data[k][j][i] + omega12->data[k][j][i - 1]) + 
                            tau * cfparams.rigth_beta / (4 * grid.hx1 * grid.hx1) * 
                            (temperature12->data[k][j][i + 1] - temperature12->data[k][j][i - 1]);

                alpha_y->data[j + 1] = B / (C - A * alpha_y->data[j]);
                beta_y->data[j + 1] = (A * beta_y->data[j] + F) / (C - A * alpha_y->data[j]);
            }

            omega->data[k + 1][grid.x2->size - 1][i] = -2 / (grid.hx2 * grid.hx2) * psi->data[k][grid.x2->size - 2][i];
            for (int j = grid.x2->size - 2; j >= 0; j--)
            {
                omega->data[k + 1][j][i] = alpha_y->data[j + 1] * omega->data[k + 1][j + 1][i] + beta_y->data[j + 1];
                omega->data[k + 1][j][grid.x1->size - 1] = -2 / (grid.hx2 * grid.hx2) * psi->data[k][j][grid.x1->size - 2];
                omega->data[k + 1][j][0] = -2 / (grid.hx2 * grid.hx2) * psi->data[k][j][1];
            }
        }

        for (int j = 1; j < grid.x2->size - 1; j++)
        {
            alpha_x->data[1] = 0;
            beta_x->data[1] = 0;

            for (int i = 1; i < grid.x1->size - 1; i++)
            {
                double A = tau * cfparams.lambda / (2 * grid.hx1 * grid.hx1);
                double B = tau * cfparams.lambda / (2 * grid.hx1 * grid.hx1);
                double C = 1 + A + B;
                double F = psi->data[k][j][i] + tau * cfparams.lambda / (2 * grid.hx2 * grid.hx2) * 
                            (psi->data[k][j + 1][i] - 2 * psi->data[k][j][i] + psi->data[k][j - 1][i]) + 
                            tau * cfparams.lambda / 2 * omega->data[k][j][i];

                alpha_x->data[i + 1] = B / (C - alpha_x->data[i] * A);
                beta_x->data[i + 1] = (F + A * beta_x->data[i]) / (C - alpha_x->data[i] * A);
            }

            psi12->data[k][j][grid.x1->size - 1] = 0; // beta_x->data[grid.x1->size - 1] / (1 - alpha_x->data[grid.x1->size - 1]);
            for (int i = grid.x1->size - 2; i >= 0; i--)
            {
                psi12->data[k][j][i] = alpha_x->data[i + 1] * psi12->data[k][j][i + 1] + beta_x->data[i + 1];
            }
        }

        for (int i = 1; i < grid.x1->size - 1; i++)
        {
            alpha_y->data[1] = 0;
            beta_y->data[1] = 0;

            for (int j = 1; j < grid.x2->size - 1; j++)
            {
                double A = tau * cfparams.lambda / (2 * grid.hx2 * grid.hx2);
                double B = tau * cfparams.lambda / (2 * grid.hx2 * grid.hx2);
                double C = 1 + A + B;
                double F = psi12->data[k][j][i] + tau * cfparams.lambda / (2 * grid.hx1 * grid.hx1) * 
                            (psi12->data[k][j][i + 1] - 2 * psi12->data[k][j][i] + psi12->data[k][j][i - 1]) + 
                            tau * cfparams.lambda / 2 * omega->data[k][j][i];

                alpha_y->data[j + 1] = B / (C - alpha_y->data[j] * A);
                beta_y->data[j + 1] = (F + A * beta_y->data[j]) / (C - alpha_y->data[j] * A);
            }

            psi->data[k + 1][grid.x2->size - 1][i] = 0; // beta_y->data[grid.x2->size - 1] / (1 - alpha_y->data[grid.x2->size - 1]);
            for (int j = grid.x2->size - 2; j >= 0; j--)
            {
                psi->data[k + 1][j][i] = alpha_y->data[j + 1] * psi->data[k + 1][j + 1][i] + beta_y->data[j + 1];
            }
        }

        for (int j = 1; j < grid.x2->size-1; j++)
        {
            for (int i = 1; i < grid.x1->size-1; i++)
            {
                u->data[j][i] = (psi->data[k][j + 1][i] - psi->data[k][j - 1][i]) / (2 * grid.hx2);
                v->data[j][i] = -(psi->data[k][j][i + 1] - psi->data[k][j][i - 1]) / (2 * grid.hx1);
            }
        }
    }
}

void helper(const char *progname)
{
    fprintf(stdout, "%s starts a calculations convective flow of fluid in 2D case.\n\n", progname);
    fprintf(stdout, "Usage:\n");
    fprintf(stdout, "    %s [OPTION]...      \n", progname);
    fprintf(stdout, "\nOptions:\n");
    fprintf(stdout, "    -G, --grashof       Reynolds number\n");
    fprintf(stdout, "    -P, --prandtl       Prandtl number\n");
    fprintf(stdout, "    -R, --reynolds      Grashof number\n");
    fprintf(stdout, "    -T, --time-points   number of nodes by time\n");
    fprintf(stdout, "    -t, --t1            \n");
    fprintf(stdout, "    -x, --x1            \n");
    fprintf(stdout, "    -X, --x-points      number of nodes by X\n");
    fprintf(stdout, "    -y, --y1            \n");
    fprintf(stdout, "    -Y, --y-points      number of nodes by Y\n");
    fprintf(stdout, "        --x-scale       \n");
    fprintf(stdout, "        --y-scale       \n");
    fprintf(stdout, "    -?, --help          \n");
}

void convective_flow_data_free()
{
    fprintf(logfile, "\nFREEE ALL ALLOCATED MEMORY\n");
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
    free_vector(grid.x1);
    free_vector(grid.x2);
    free_vector(t);
}

void save_data()
{
    fprintf(logfile, "SAVE DATA\n");
    FILE *temperature_file = fopen("./output_data/temperature.txt", "w");
    FILE *omega_file = fopen("./output_data/omega.txt", "w");
    FILE *psi_file = fopen("./output_data/psi.txt", "w");
    FILE *velocity_u = fopen("./output_data/u_vel.txt", "w");
    FILE *velocity_v = fopen("./output_data/v_vel.txt", "w");
    if (temperature_file == NULL || omega_file == NULL || psi_file == NULL)
        return;

    for (int j = 0; j < grid.x2->size; j++)
    {
        for (int i = 0; i < grid.x1->size; i++)
        {
            fprintf(temperature_file, "%.4f ", temperature->data[K - 1][j][i]);
            fprintf(omega_file, "%.4f ", omega->data[K - 1][j][i]);
            fprintf(psi_file, "%.4f ", psi->data[K - 1][j][i]);
            fprintf(velocity_u, "%.4f ", u->data[j][i]);
            fprintf(velocity_v, "%.4f ", v->data[j][i]);
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
    fprintf(logfile, "\tSuccessful!\n");

    fprintf(logfile, "SAVE TIME DATA\n");
    FILE *temperature_time_data = fopen("./output_data/temperature_time_data.txt", "w");
    if (temperature_time_data == NULL)
        return;

    for (int k = 0; k < K; k += 10)
    {
        for (int j = 0; j < grid.x2->size; j++)
            for (int i = 0; i < grid.x1->size; i++)
                fprintf(temperature_time_data, "%.4f ", temperature->data[k][j][i]);
        fprintf(temperature_time_data, "%c", '\n');
    }
    fclose(temperature_time_data);
}

double getMaxAbs(Matrix3D *p)
{
    double max = 0.0;
    double maxitem = 0.0;
    double tmp = 0.0;
    for (int i = 0; i < grid.x1->size; i++)
    {
        for (int j = 0; j < grid.x2->size; j++)
        {
            tmp = fabs(p->data[K - 1][j][i] - p->data[K - 2][j][i]);
            if (max <= tmp)
                max = tmp;

            if (maxitem <= fabs(p->data[K - 1][j][i]))
                maxitem = p->data[K - 1][j][i];
        }
    }

    fprintf(logfile, "MAX PSI ITEM: %f\n", maxitem);

    return max;
}

void inititalize_convective_flow_environment()
{
    t = init_vector(K);
    grid.x1 = init_vector(N);
    grid.x2 = init_vector(M);
    alpha_x = init_vector(N);
    alpha_y = init_vector(M);
    beta_x = init_vector(N);
    beta_y = init_vector(M);

    u = init_matrix(N, M);
    v = init_matrix(N, M);

    psi = init_matrix3d(N, M, K);
    psi12 = init_matrix3d(N, M, K);
    omega = init_matrix3d(N, M, K);
    omega12 = init_matrix3d(N, M, K);
    temperature = init_matrix3d(N, M, K);
    temperature12 = init_matrix3d(N, M, K);
}

int main(int argc, char *argv[])
{
    logfile = fopen("logfile", "a");
    double x_1 = 1.0;
    double y_1 = 1.0;
    double t_1 = 1.0;

    if (argc <= 1)
    {
        fprintf(logfile, "No arguments are given!\n\n");
        helper(argv[0]);
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
        {NULL, 0, NULL, 0}};
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
                fprintf(logfile, "Wrong number of X steps: %d\n", N);
                exit(1);
            }
            break;
        case 'Y':
            M = atoi(optarg);
            if (M < 10)
            {
                fprintf(logfile, "Wrong number of Y steps: %d\n", M);
                exit(1);
            }
            break;
        case 'T':
            K = atoi(optarg);
            if (K < 10)
            {
                fprintf(logfile, "Wrong number of time steps: %d\n", K);
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
            fprintf(stdout, "Try \"./convective_flow --help\" for more information.\n");
            exit(0);
        }
    }

    double hx;
    double hy;

    K = K * t_1 + 1;
    N = N * x_1 + 1;
    M = M * y_1 + 1;

    inititalize_convective_flow_environment();

    if (alloc_error)
    {
        fprintf(logfile, "Failed to allocate memory\n");
        convective_flow_data_free();
        exit(1);
    }

    tau = t_1 / (K - 1);
    grid.hx1 = x_1 / (grid.x1->size - 1);
    grid.hx2 = y_1 / (grid.x2->size - 1);
    hx = grid.hx1;
    hy = grid.hx2;

    cfparams.xi = 1.0 / (cfparams.Re * cfparams.Pr);
    cfparams.nu = 1.0 / cfparams.Re;
    cfparams.rigth_beta = cfparams.Gr / (cfparams.Re * cfparams.Re);
    cfparams.lambda = 0.5;

    for (int i = 0; i < grid.x1->size; i++)
        grid.x1->data[i] = i * hx;

    for (int j = 0; j < grid.x2->size; j++)
        grid.x2->data[j] = j * hy;

    for (int k = 0; k < K; k++)
        t->data[k] = k * tau;

    for (int j = 0; j < grid.x2->size; j++)
    {
        for (int i = 0; i < grid.x1->size; i++)
        {
            temperature->data[0][j][i] = T0;
            temperature12->data[0][j][i] = T0;
        }

        temperature->data[0][j][grid.x1->size-1] = T1;
        temperature12->data[0][j][grid.x1->size-1] = T1;
    }

    fprintf(logfile, "RUN CALCULATIONS\n");
    calculations();
    fprintf(logfile, "\tSuccessful!\n");

    save_data();

    fprintf(logfile, "%f\n", getMaxAbs(psi));
    fprintf(logfile, "\tSuccessful!\n");
    convective_flow_data_free();
    fprintf(logfile, "\tSuccessful!\n");

    return 0;
}
