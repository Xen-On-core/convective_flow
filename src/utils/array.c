#include <malloc.h>
#include "utils/array.h"

float *init_zero_1d(int size) {
    float *array = malloc(sizeof(float) * size);
    if (array == NULL) {
        return NULL;
    } 
    for (int i = 0; i < size; i++)
        array[i] = 0;

    return array;
}

float **init_zero_2d(int x_size, int y_size) {
    float **array = malloc(sizeof(float*) * y_size);
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

float ***init_zero_3d(int x_size, int y_size, int z_size) {
    float ***array = malloc(sizeof(float**) * z_size);
    if (array == NULL) {
        return NULL;
    }
    
    for (int k = 0; k < z_size; k++)
    {
        array[k] = malloc(sizeof(float*) * y_size);
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

void free_1d_array(float *array) {
    if (array != NULL) {
        free(array);
    }
}

void free_2d_array(float **array, int M) {
    if (array != NULL) {
        for (int j = 0; j < M; j++)
            if (array[j] != NULL) {
                free(array[j]);
            }
        free(array);
    }
}

void free_3d_array(float ***array, int K, int M) {
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

void print_array_1d(float *array, int size) {
    for (int i = 0; i < size; i++)
        printf("%.4f ", array[i]);
    printf("\n");
}

void print_array_2d(float **array, int M, int N) {
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++)
            printf("%.4f ", array[j][i]);
        printf("\n");
    }
    printf("\n");
}

void print_array_3d(float ***array, int K, int M, int N) {
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