#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

int N;
int M;
int K;

Vector init_zero_vector(int size) {
    Vector array = malloc(sizeof(float) * size);
    if (array == NULL) {
        return NULL;
    } 
    for (int i = 0; i < size; i++)
        array[i] = 0;

    return array;
}

Matrix init_zero_matrix(int x_size, int y_size) {
    Matrix array = malloc(sizeof(float*) * y_size);
    if (array == NULL) {
        return NULL;
    }
    
    for (int j = 0; j < y_size; j++){
        array[j] = malloc(sizeof(float) * x_size);
        if (array[j] == NULL) {
            return NULL;
        }
    }
    
    for (int j = 0; j < y_size; j++)
        for (int i = 0; i < x_size; i++)
            array[j][i] = 0.0;

    return array;
}

Array3D init_zero_array3d(int x_size, int y_size, int z_size) {
    Array3D array = malloc(sizeof(Matrix) * z_size);
    if (array == NULL) {
        return NULL;
    }
    
    for (int k = 0; k < z_size; k++)
    {
        array[k] = malloc(sizeof(Vector) * y_size);
        if (array[k] == NULL) {
            return NULL;
        }
        for (int j = 0; j < y_size; j++){
            array[k][j] = malloc(sizeof(float) * x_size);
            if (array[k][j] == NULL) {
                return NULL;
            }
        }
    }


    for (int k = 0; k < z_size; k++)
        for (int j = 0; j < y_size; j++)
            for (int i = 0; i < x_size; i++)
                array[k][j][i] = 0.0;

    return array;
}

void free_vector(Vector array) {
    if (array != NULL) {
        free(array);
    }
}

void free_matrix(Matrix array) {
    if (array != NULL) {
        for (int j = 0; j < M; j++)
            if (array[j] != NULL) {
                free(array[j]);
            }
        free(array);
    }
}

void free_array3d(Array3D array) {
    if (array != NULL) {
        for (int k = 0; k < K; k++){
            for (int j = 0; j < M; j++){
                if (array[k][j] != NULL) {
                    free(array[k][j]);
                }
            }
        }
        free(array);
    }
}

void print_vector(float *array, int size) {
    for (int i = 0; i < size; i++)
        printf("%.4f ", array[i]);
    printf("\n");
}

void print_matrix(float **array) {
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++)
            printf("%.4f ", array[j][i]);
        printf("\n");
    }
    printf("\n");
}

void print_array3d(float ***array) {
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