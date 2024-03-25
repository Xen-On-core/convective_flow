#include <stdio.h>
#include <stdlib.h>
#include "utils/array.h"

extern int N;
extern int M;
extern int K;
extern int T;

Vector init_zero_vector(int size)
{
    Vector vector = malloc(sizeof(double) * size);
    if (vector == NULL)
        return NULL;
 
    for (int i = 0; i < size; i++)
        vector[i] = 0;

    return vector;
}

Matrix init_zero_matrix(int x_size, int y_size)
{
    Matrix matrix = malloc(sizeof(Vector) * y_size);
    if (matrix == NULL)
        return NULL;

    for (int j = 0; j < y_size; j++)
        matrix[j] = init_zero_vector(x_size);

    return matrix;
}

Matrix3D init_zero_matrix3d(int x_size, int y_size, int z_size)
{
    Matrix3D matrix = malloc(sizeof(Matrix) * z_size);
    if (matrix == NULL)
        return NULL;
    
    for (int k = 0; k < z_size; k++)
        matrix[k] = init_zero_matrix(x_size, y_size);

    return matrix;
}

TimeVector init_zero_time_vector(int time_size, int size)
{
    TimeVector vector = malloc(sizeof(Vector) * time_size);
    if (vector == NULL)
        return NULL;
 
    for (int t = 0; t < time_size; t++)
        vector[t] = init_zero_vector(size);

    return vector;
}

TimeMatrix init_zero_time_matrix(int time_size, int x_size, int y_size)
{
    TimeMatrix matrix = malloc(sizeof(Matrix) * time_size);
    if (matrix == NULL)
        return NULL;

    for (int t = 0; t < time_size; t++)
        matrix[t] = init_zero_matrix(x_size, y_size);

    return matrix;
}

TimeMatrix3D init_zero_time_matrix3d(int time_size, int x_size, int y_size, int z_size)
{
    TimeMatrix3D matrix = malloc(sizeof(Matrix3D) * time_size);
    if (matrix == NULL)
        return NULL;
    
    for (int t = 0; t < time_size; t++)
        matrix[t] = init_zero_matrix3d(x_size, y_size, z_size);

    return matrix;
}

void free_vector(Vector vector)
{
    if (vector != NULL)
        free(vector);
}

void free_matrix(Matrix matrix)
{
    if (matrix == NULL)
        return;
    
    for (int j = 0; j < M; j++)
        free_vector(matrix[j]);
    free(matrix);
}

void free_matrix3d(Matrix3D matrix)
{
    if (matrix == NULL)
        return;
    
    for (int k = 0; k < K; k++)
        free_matrix(matrix[k]);
    free(matrix);
}

void free_time_vector(TimeVector vector)
{
    if (vector == NULL)
        return;

    for (int t = 0; t < T; t++)
        free_vector(vector[t]);
    free(vector);
}

void free_time_matrix(TimeMatrix matrix)
{
    if (matrix == NULL)
        return;

    for (int t = 0; t < T; t++)
        free_matrix(matrix[t]);
    free(matrix);
}

void free_time_matrix3d(TimeMatrix3D matrix)
{
    if (matrix == NULL)
        return;

    for (int t = 0; t < T; t++)
        free_matrix3d(matrix[t]);   
    free(matrix);
}

void print_vector(Vector vector)
{
    for (int i = 0; i < N; i++)
        printf("%.4f ", vector[i]);
    printf("\n");
}

void print_matrix(Matrix matrix) {
    for (int j = 0; j < M; j++)
            print_vector(matrix[j]);
    printf("\n");
}

void print_matrix3d(Matrix3D matrix) {
    for (int k = 0; k < K; k++)
        print_matrix(matrix[k]);
    printf("\n");
}

void print_time_vector(TimeVector vector, int step) {
    if (step == 0)
    {
        printf("ERROR! Step for time cannot be equal zero!");
        return;
    }

    for(int t = 0; t < T; t+=step)
        print_vector(vector[t]);
    printf("\n");
}

void print_time_matrix(TimeMatrix matrix, int step) {
    if (step == 0)
    {
        printf("ERROR! Step for time cannot be equal zero!");
        return;
    }

    for (int t = 0; t < T; t+=step)
        print_matrix(matrix[t]);
    printf("\n");
}

void print_time_matrix3d(TimeMatrix3D matrix, int step) {
    if (step == 0)
    {
        printf("ERROR! Step for time cannot be equal zero!");
        return;
    }

    for (int t = 0; t < T; t+=step) 
        print_matrix3d(matrix[t]);
    printf("\n");
}