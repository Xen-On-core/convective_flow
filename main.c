#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Инициализируем следующие переменные, описывающие сетку для вычислений:
 *      N : кол-во узлов по X
 *      M : кол-во узлов по Y
 *      K : кол-во узлов по времени
 *      x0, y0 : значение по X и Y на удаленных от осей координат границах
 *      t0 : гранчное время
 *      hx, hy, tau : шаг по X, Y и времени соответственно
 *      
 */
#define N ( 40 )
#define M ( 40 )
#define K ( 10000 )
#define x0 ( 1.0 )
#define y0 ( 1.0 )
#define t0 ( 1.0 )
#define hx ( x0 / N )
#define hy ( y0 / M )
#define tau ( t0 / K )

/*
 * Инициализируем константы:
 *      Re : число Рейнольдса
 *      Pr : число Прандтля
 *      Gr : число Грасгофа
 *      xi : коэфф. температуропроводности
 *      nu : коэфф. кинематической вязкости
 *      g : ускорение силы тяжести
 *      right_beta : коэфф. температурного расширения
 */
#define Re ( 1.0 )
#define Gr ( 10000.0 )
#define Pr ( 1.0 )
#define xi ( 1 / (Re * Pr) )
#define nu ( 1 / Re )
#define g ( 9.8 )
#define rigth_beta ( Gr / Re)

float x[N];
float y[M];
float t[K];

float T[K][M][N];
float T12[K][M][N];
float omega[K][M][N];
float omega12[K][M][N];
float psi[M][N];
float u[M][N];
float v[M][N];

float alpha_x[N];
float beta_x[N];
float alpha_y[M];
float beta_y[M];

void init_zero_1d(float array[], int size) {
    for (int i = 0; i < size; i++)
        array[i] = 0;
}

void init_zero_2d(float array[][N]) {
    for (int j = 0; j < M; j++)
        for (int i = 0; i < N; i++)
            array[j][i] = 0;
}

void init_zero_3d(float array[][M][N]) {
    for (int k = 0; k < K; k++)
        for (int j = 0; j < M; j++)
            for (int i = 0; i < N; i++)
                array[k][j][i] = 0;
}

void print_array_1d(float array[], int size) {
    for (int i = 0; i < size; i++)
        printf("%.4f ", array[i]);
    printf("\n");
}

void print_array_2d(float array[][N]) {
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++)
            printf("%.4f ", array[j][i]);
        printf("\n");
    }
    printf("\n");
}

void print_array_3d(float array[][M][N]) {
    for (int k = 0; k < K; k++) {
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < N; i++)
                printf("%.4f ", array[k][j][i]);
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

void init_boundary_conditions(float array[K][M][N], char side, int a, int b, float value) {
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
    for (int k = 0; k < K; k++)
    {
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

        /*
         * Run back-for by X for over Y on n+1/2 layer
         */
        float coef_x = tau/(2*hx);
        for (int j = 1; j < M-1; j++)
        {
            alpha_x[1] = 1;
            beta_x[1] = 0;

            for (int i = 1; i < N-1; i++)
            {
                float A = coef_x * (xi/hx + u[j][i]/2);
                float B = coef_x * (xi/hx - u[j][i]/2);
                float C = 1 + tau * xi / (hx * hx);
                float F = T[k][j][i] - tau*v[j][i]/(4*hy)*(T[k][j+1][i] - T[k][j-1][i]) + \
                            tau*xi/(2*hy*hy)*(T[k][j+1][i] - 2*T[k][j][i] + T[k][j-1][i]);

                alpha_x[i+1] = B / (C - A * alpha_x[i]);
                beta_x[i+1] = (A * beta_x[i] + F) / (C - A * alpha_x[i]);
            }

            T12[k][j][N-1] = beta_x[N-1]/(1-alpha_x[N-1]);
            for (int i = N-2; i >= 0; i--)
            {
                T12[k][j][i] = alpha_x[i+1] * T12[k][j][i+1] + beta_x[i+1];
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
            alpha_y[1] = 0;
            beta_y[1] = 1.0/8.0;

            for (int j = 1; j < M-1; j++)
            {
                float A = coef_y * (xi/hy + v[j][i]/2);
                float B = coef_y * (xi/hy - v[j][i]/2);
                float C = 1 + tau * xi / (hy * hy);
                float F = T12[k][j][i] - tau*u[j][i]/(4*hx)*(T12[k][j][i+1] - T12[k][j][i-1]) + \
                            tau*xi/(2*hx*hx)*(T12[k][j][i+1] - 2*T12[k][j][i] + T12[k][j][i-1]);

                alpha_y[j+1] = B / (C - A * alpha_y[j]);
                beta_y[j+1] = (A * beta_y[j] + F) / (C - A * alpha_y[j]);
            }

            T[k+1][M-1][i] = 1.0;
            for (int j = M-2; j >= 0; j--)
            {
                T[k+1][j][i] = alpha_y[j+1] * T[k+1][j+1][i] + beta_y[j+1];
                T[k+1][j][0] = T[k+1][j+1][1];
                T[k+1][j][N-1] = T[k+1][j+1][N-2];

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
                            tau*rigth_beta/(4*hx*hx)*(T[k][j][i+1] - T[k][j][i-1]);

                alpha_x[i+1] = B / (C - A * alpha_x[i]);
                beta_x[i+1] = (A * beta_x[i] + F) / (C - A * alpha_x[i]);
            }

            omega12[k][j][N-1] = -2 / (hx*hx) * psi[j][N-2];
            for (int i = N-2; i >= 0; i--)
                omega12[k][j][i] = alpha_x[i+1] * omega12[k][j][i+1] + beta_x[i+1];
            omega12[k][j][0] = -2 / (hx*hx) * psi[j][1];
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
                            tau*rigth_beta/(4*hx*hx)*(T[k][j][i+1] - T[k][j][i-1]);

                alpha_y[j+1] = B / (C - A * alpha_y[j]);
                beta_y[j+1] = (A * beta_y[j] + F) / (C - A * alpha_y[j]);
            }

            omega[k+1][M-1][i] = -2 / (hy*hy) * psi[M-2][i];
            for (int j = M-2; j >= 0; j--)
                omega[k+1][j][i] = alpha_y[j+1] * omega[k+1][j+1][i] + beta_y[j+1];
            omega[k+1][0][i] = -2 / (hy*hy) * psi[1][i];
        }

        /*
         * Calculations of PSI with star-scheme
         */
        float psi_coef = (hx*hx*hy*hy/(2*hx*hx + 2*hy*hy));
        for (int j = 1; j < M-1; j++)
            for (int i = 1; i < N-1; i++)
            {
                psi[j][i] = psi_coef * (
                                (psi[j][i+1] + psi[j][i-1])/(hx*hx) + \
                                (psi[j+1][i] + psi[j-1][i])/(hy*hy) + \
                                omega[k][j][i]
                            );
            }
    }
}

int main(int argc, char *argv[]) {
    
    init_zero_1d(x, N);
    init_zero_1d(y, M);
    init_zero_1d(t, K);

    for (int i = 0; i < N; i++)
        x[i] = i*hx;

    for (int j = 0; j < M; j++)
        y[j] = j*hy;

    for (int k = 0; k < K; k++)
        t[k] = k*tau;

    init_zero_2d(u);
    init_zero_2d(v);
    init_zero_3d(T);
    init_zero_3d(T12);
    init_zero_3d(omega);
    init_zero_3d(omega12);
    init_zero_2d(psi);
    init_zero_1d(alpha_x, N);
    init_zero_1d(beta_x, N);
    init_zero_1d(alpha_y, M);
    init_zero_1d(beta_y, M);

    init_boundary_conditions(T, 'b', 0, 1, 1.0/8.0);
    init_boundary_conditions(T, 't', 0, 1, 1.0);
    init_boundary_conditions(T12, 'b', 0, 1, 1.0/8.0);
    init_boundary_conditions(T12, 't', 0, 1, 1.0);

    calculations();
    
    FILE* temp = fopen("./Temp.txt", "w");
    FILE* Omg = fopen("./Omg.txt", "w");
    FILE* Psi = fopen("./Psi.txt", "w");
    if (temp == NULL || Omg == NULL || Psi == NULL)
        return -1;

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            fprintf(temp, "%.4f ", T[K-1][j][i]);
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

    FILE* temp_all = fopen("./temp_all.txt", "w");
    if (temp_all == NULL)
        return -1;

    for (int k = 0; k < K; k+=10)
    {
        for (int j = 0; j < M; j++)
        {
            for (int i = 0; i < N; i++)
            {
                fprintf(temp_all, "%.4f ", T[k][j][i]);
            }
        }
        fprintf(temp_all, "%c", '\n');
    }

    fclose(temp_all);

    // print_array_2d(T12[K-1]);
    // print_array_2d(T[K-1]);
    // print_array_2d(omega[K-1]);
    // print_array_2d(psi);
    // print_array_2d(u);
    // print_array_2d(v);

    // print_array_3d(T);
    // print_array_3d(T12);
    // print_array_3d(omega);
    // print_array_3d(omega12);
    // print_array_2d(u);
    // print_array_2d(v);

    return 0;
}