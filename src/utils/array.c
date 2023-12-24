#include <stdio.h>
#include <stdlib.h>
#include "utils/array.h"

int N;
int M;
int K;

Vector init_zero_vector(int size) {
    Vector vector = malloc(sizeof(float) * size);
    if (vector == NULL) {
        return NULL;
    } 
    for (int i = 0; i < size; i++)
        vector[i] = 0;

    return vector;
}

Matrix init_zero_matrix(int x_size, int y_size) {
    Matrix matrix = malloc(sizeof(float*) * y_size);
    if (matrix == NULL) {
        return NULL;
    }
    
    for (int j = 0; j < y_size; j++){
        matrix[j] = malloc(sizeof(float) * x_size);
        if (matrix[j] == NULL) {
            return NULL;
        }
    }
    
    for (int j = 0; j < y_size; j++)
        for (int i = 0; i < x_size; i++)
            matrix[j][i] = 0.0;

    return matrix;
}

Matrix3D init_zero_matrix3d(int x_size, int y_size, int z_size) {
    Matrix3D matrix = malloc(sizeof(Matrix) * z_size);
    if (matrix == NULL) {
        return NULL;
    }
    
    for (int k = 0; k < z_size; k++)
    {
        matrix[k] = malloc(sizeof(Vector) * y_size);
        if (matrix[k] == NULL) {
            return NULL;
        }
        for (int j = 0; j < y_size; j++){
            matrix[k][j] = malloc(sizeof(float) * x_size);
            if (matrix[k][j] == NULL) {
                return NULL;
            }
        }
    }


    for (int k = 0; k < z_size; k++)
        for (int j = 0; j < y_size; j++)
            for (int i = 0; i < x_size; i++)
                matrix[k][j][i] = 0.0;

    return matrix;
}

void free_vector(Vector vector) {
    if (vector != NULL) {
        free(vector);
    }
}

void free_matrix(Matrix matrix) {
    if (matrix != NULL) {
        for (int j = 0; j < M; j++)
            if (matrix[j] != NULL) {
                free(matrix[j]);
            }
        free(matrix);
    }
}

void free_matrix3d(Matrix3D matrix) {
    if (matrix != NULL) {
        for (int k = 0; k < K; k++){
            for (int j = 0; j < M; j++){
                if (matrix[k][j] != NULL) {
                    free(matrix[k][j]);
                }
            }
        }
        free(matrix);
    }
}

void print_vector(Vector vector, int size) {
    for (int i = 0; i < size; i++)
        printf("%.4f ", vector[i]);
    printf("\n");
}

void print_matrix(Matrix matrix) {
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++)
            printf("%.4f ", matrix[j][i]);
        printf("\n");
    }
    printf("\n");
}

void print_matrix3d(Matrix3D matrix) {
    for (int k = 0; k < K; k++) {
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < N; i++)
                printf("%.4f ", matrix[k][j][i]);
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}