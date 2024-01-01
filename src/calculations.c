#include "structures.h"
#include "utils/array.h"
#include <malloc.h>

int NewBounds(Bounds *B, float t_1, float x_1, float y_1) {
    B->t_1 = 1.0;
    B->x_1 = 1.0;
    B->y_1 = 1.0;

    B->x_0 = 0.0;
    B->y_0 = 0.0;
    return 0;
}

int NewSteps(Steps *S, Bounds *B, int N, int M, int K) { 
    S->K = K;
    S->M = M;
    S->N = N;
    S->hx = (B->x_1) / ((float)N);
    S->hy = (B->y_1) / ((float)M);
    S->tau = (B->t_1) / ((float)K);
    return 0;
}

int NewConstants(Constants *C) {
    C->Gr = 10000.0;
    C->Re = 1.0;
    C->Pr = 1.0;
    C->xi = 1.0 / (C->Re * C->Pr);
    C->nu = 1.0 / (C->Re);
    C->rigth_beta = C->Gr / (C->Re * C->Re);
    return 0;
}

int AllocSolutionArrays(Solution *SLTB, Steps *STP) {
    SLTB->x = init_row(STP->N+1);
    SLTB->y = init_row(STP->M+1);
    SLTB->t = init_row(STP->K);
    SLTB->temperature = init_tensor(STP->N, STP->M, STP->K);
    SLTB->temperature12 = init_tensor(STP->N, STP->M, STP->K);
    SLTB->omega = init_tensor(STP->N, STP->M, STP->K);
    SLTB->omega12 = init_tensor(STP->N, STP->M, STP->K);
    SLTB->psi = init_matrix(STP->N, STP->M);
    SLTB->u = init_matrix(STP->N, STP->M);
    SLTB->v = init_matrix(STP->N, STP->M);
    SLTB->alpha_x = init_row(STP->N);
    SLTB->beta_x = init_row(STP->N);
    SLTB->alpha_y = init_row(STP->M);
    SLTB->beta_y = init_row(STP->M);
    if ( SLTB->x == NULL || SLTB->y == NULL || SLTB->t == NULL || \
        SLTB->temperature == NULL || SLTB->temperature12 == NULL || \
        SLTB->omega == NULL || SLTB->omega12 == NULL || SLTB->psi == NULL || \
        SLTB->u == NULL || SLTB->v == NULL || SLTB->alpha_x == NULL || \
        SLTB->beta_y == NULL )
    {
        printf("Failed to allocate memory\n");
        FreeSLTButionContainers(SLTB, STP);
        return 1;
    }

    // TODO: Replace it in right place
    for (int i = 0; i < STP->N+1; i++)
        SLTB->x[i] = i*STP->hx;

    for (int j = 0; j < STP->M+1; j++)
        SLTB->y[j] = j*STP->hy;

    for (int k = 0; k < STP->K; k++)
        SLTB->t[k] = k*STP->tau;

    


    return 0;
}



void init_boundary_conditions(Matrix3D array, char side, int a, int b, float value, Steps* STP) {
    int N = STP->N, M = STP->M, K = STP->K ;

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


void calculations(Solution *SLTB, double T_0, double T_1) {
    int N = SLTB->STP->N, M = SLTB->STP->M, K = SLTB->STP->K;
    float hx = SLTB->STP->hx, hy = SLTB->STP->hy, tau = SLTB->STP->tau; 
    float xi = SLTB->Constants->xi,
          Gr = SLTB->Constants->Gr,
          Pr = SLTB->Constants->Pr,
          nu = SLTB->Constants->nu, 
          rigth_beta = SLTB->Constants->rigth_beta;



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
            SLTB->alpha_x[1] = 0;
            SLTB->beta_x[1] = T_0;

            for (int i = 1; i < N-1; i++)
            {
                float A = coef_x * (xi/hx + SLTB->u[j][i]/2);
                float B = coef_x * (xi/hx - SLTB->u[j][i]/2);
                float C = 1 + tau * xi / (hx * hx);
                float F = SLTB->temperature[k][j][i] - tau*SLTB->v[j][i]/(4*hy)*(SLTB->temperature[k][j+1][i] - SLTB->temperature[k][j-1][i]) + \
                            tau*xi/(2*hy*hy)*(SLTB->temperature[k][j+1][i] - 2*SLTB->temperature[k][j][i] + SLTB->temperature[k][j-1][i]);
                

                SLTB->alpha_x[i+1] = B / (C - A * SLTB->alpha_x[i]);
                SLTB->beta_x[i+1] = (A * SLTB->beta_x[i] + F) / (C - A * SLTB->alpha_x[i]);
            }

            SLTB->temperature12[k][j][N-1] = T_1;
            for (int i = N-2; i >= 0; i--)
            {
                SLTB->temperature12[k][j][i] = SLTB->alpha_x[i+1] * SLTB->temperature12[k][j][i+1] + SLTB->beta_x[i+1];
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
            SLTB->alpha_y[1] = 1;
            SLTB->beta_y[1] = 0;

            for (int j = 1; j < M-1; j++)
            {
                float A = coef_y * (xi/hy + SLTB->v[j][i]/2);
                float B = coef_y * (xi/hy - SLTB->v[j][i]/2);
                float C = 1 + tau * xi / (hy * hy);
                float F = SLTB->temperature12[k][j][i] - tau*SLTB->u[j][i]/(4*hx)*(SLTB->temperature12[k][j][i+1] - SLTB->temperature12[k][j][i-1]) + \
                            tau*xi/(2*hx*hx)*(SLTB->temperature12[k][j][i+1] - 2*SLTB->temperature12[k][j][i] + SLTB->temperature12[k][j][i-1]);

                SLTB->alpha_y[j+1] = B / (C - A * SLTB->alpha_y[j]);
                SLTB->beta_y[j+1] = (A * SLTB->beta_y[j] + F) / (C - A * SLTB->alpha_y[j]);
            }

            SLTB->temperature[k+1][M-1][i] = SLTB->beta_y[N-1]/(1-SLTB->alpha_y[N-1]);
            for (int j = M-2; j >= 0; j--)
            {
                SLTB->temperature[k+1][j][i] = SLTB->alpha_y[j+1] * SLTB->temperature[k+1][j+1][i] + SLTB->beta_y[j+1];
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
            SLTB->alpha_x[1] = 0;
            SLTB->beta_x[1] = -2 / (hx*hx) * SLTB->psi[j][1];

            for (int i = 1; i < N-1; i++)
            {
                float A = coef_x * (nu/hx + SLTB->u[j][i]/2);
                float B = coef_x * (nu/hx - SLTB->u[j][i]/2);
                float C = 1 + tau * nu / (hx * hx);
                float F = omega[k][j][i] - tau*SLTB->v[j][i]/(4*hy)*(omega[k][j+1][i] - omega[k][j-1][i]) + \
                            tau*nu/(2*hy*hy)*(omega[k][j+1][i] - 2*omega[k][j][i] + omega[k][j-1][i]) + \
                            tau*rigth_beta/(4*hx*hx)*(SLTB->temperature[k][j][i+1] - SLTB->temperature[k][j][i-1]);

                SLTB->alpha_x[i+1] = B / (C - A * SLTB->alpha_x[i]);
                SLTB->beta_x[i+1] = (A * SLTB->beta_x[i] + F) / (C - A * SLTB->alpha_x[i]);
            }

            SLTB->omega12[k][j][N-1] = -2 / (hx*hx) * SLTB->psi[j][N-2];
            for (int i = N-2; i >= 0; i--)
                SLTB->omega12[k][j][i] = SLTB->alpha_x[i+1] * SLTB->omega12[k][j][i+1] + SLTB->beta_x[i+1];
        }

        /*
         * Run back-for by Y for over X on n+1 layer
         */
        for (int i = 1; i < N-1; i++)
        {
            // alpha - a second-order condition.
            // beta - a first-order condition.
            SLTB->alpha_y[1] = 0;
            SLTB->beta_y[1] = -2 / (hy*hy) * SLTB->psi[1][i];

            for (int j = 1; j < M-1; j++)
            {
                float A = coef_y * (nu/hy - SLTB->v[j][i]/2);
                float B = coef_y * (nu/hy + SLTB->v[j][i]/2);
                float C = 1 + tau * nu / (hy * hy);
                float F = SLTB->omega12[k][j][i] - tau*SLTB->u[j][i]/(4*hx)*(SLTB->omega12[k][j][i+1] - SLTB->omega12[k][j][i-1]) + \
                            tau*nu/(2*hx*hx)*(SLTB->omega12[k][j][i+1] - 2*SLTB->omega12[k][j][i] + SLTB->omega12[k][j][i-1]) + \
                            tau*rigth_beta/(4*hx*hx)*(SLTB->temperature[k][j][i+1] - SLTB->temperature[k][j][i-1]);

                SLTB->alpha_y[j+1] = B / (C - A * SLTB->alpha_y[j]);
                SLTB->beta_y[j+1] = (A * SLTB->beta_y[j] + F) / (C - A * SLTB->alpha_y[j]);
            }

            SLTB->omega[k+1][M-1][i] = -2 / (hy*hy) * SLTB->psi[M-2][i];
            for (int j = M-2; j >= 0; j--)
            {
                SLTB->omega[k+1][j][i] = SLTB->alpha_y[j+1] * omega[k+1][j+1][i] + SLTB->beta_y[j+1];
                SLTB->omega[k+1][j][N-1] = -2 / (hy*hy) * SLTB->psi[j][N-2];
                SLTB->omega[k+1][j][0] = -2 / (hy*hy) * SLTB->psi[j][1];
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
                SLTB->psi[j][i] = psi_coef * (
                                (SLTB->psi[j][i+1] + SLTB->psi[j][i-1])/(hx*hx) + \
                                (SLTB->psi[j+1][i] + SLTB->psi[j-1][i])/(hy*hy) + \
                                SLTB->omega[k][j][i]);
            }
        }

        /*
         * Calculations of velocities
         */
        for (int j = 1; j < M-1; j++)
        {
            for (int i = 1; i < N-1; i++)
            {
                SLTB->u[j][i] = (SLTB->psi[j+1][i] - SLTB->psi[j-1][i]) / (2 * hy);
                SLTB->v[j][i] = -(SLTB->psi[j][i+1] - SLTB->psi[j][i-1]) / (2 * hx);
            }
        }
    }
}


int NewSolution(StarterPack STRPCK){
    int exitcode = 0;
    Bounds* BNDS = (Bounds*)malloc(sizeof(Bounds));
    Steps* STP = (Steps*)malloc(sizeof(Steps));
    Constants *CNST = (Constants *)malloc(sizeof(Constants));
    exitcode = NewBounds(BNDS, STRPCK.t_1, STRPCK.x_1, STRPCK.y_1);
    if ( exitcode != 0 ){
        return exitcode;
    }
    exitcode = NewSteps(STP, BNDS,STRPCK.N,STRPCK.M,STRPCK.K);
    if ( exitcode != 0 ){
        return exitcode;
    }
    exitcode = NewConstants(CNST);
    if ( exitcode != 0 ){
        return exitcode;
    }
    Solution * SLTB = \
        (Solution *)malloc(sizeof(Solution));

    exitcode = AllocSolutionArrays(SLTB, STP);
    if ( exitcode != 0 ){
        return exitcode;
    }
    SLTB->STP = STP;
    SLTB->Bounds = BNDS;
    SLTB->Constants = CNST;
    
    for (int j = 0; j < STP->M; j++)
    {
        for (int i = 0; i < STP->N; i++)
        {
            SLTB->temperature[0][j][i] = exp(-190*(SLTB->x[i])) + STRPCK.T_0;
            SLTB->temperature12[0][j][i] = exp(-190*(SLTB->x[i])) + STRPCK.T_0;
        }
    }
    init_boundary_conditions(SLTB->temperature, 'l', 0, 1, STRPCK.T_0, STP);
    init_boundary_conditions(SLTB->temperature, 'r', 0, 1, STRPCK.T_1, STP);
    init_boundary_conditions(SLTB->temperature12, 'l', 0, 1, STRPCK.T_0, STP);
    init_boundary_conditions(SLTB->temperature12, 'r', 0, 1, STRPCK.T_0, STP);



}