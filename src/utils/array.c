#include <stdio.h>
#include <stdlib.h>
#include "utils/array.h"

_vector
init_zero_vector(int size)
{
    _vector data = malloc(sizeof(double) * size);
    if (data == NULL)
        return NULL;

    for (int i = 0; i < size; i++)
        data[i] = 0;
    
    return data;
}

_matrix
init_zero_matrix(int x_size, int y_size)
{
    _matrix data = malloc(sizeof(_matrix) * y_size);
    if (data == NULL)
        return NULL;

    for (int i = 0; i < y_size; i++)
        data[i] = init_zero_vector(x_size);

    return data;
}

_matrix3D
init_zero_matrix3d(int x_size, int y_size, int z_size)
{
    _matrix3D data = malloc(sizeof(_matrix3D) * z_size);
    if (data == NULL)
        return NULL;
    
    for (int i = 0; i < z_size; i++)
        data[i] = init_zero_matrix(x_size, y_size);

    return data;
}

Vector *
init_vector(int size)
{
    Vector *vector = (Vector *) malloc(sizeof(Vector));

    vector->size = size;
    vector->data = init_zero_vector(size);

    return vector;
}

Matrix *
init_matrix(int x_size, int y_size)
{
    Matrix *matrix = (Matrix *) malloc(sizeof(Matrix));
    
    matrix->y_size = y_size;
    matrix->x_size = x_size;
    matrix->data = init_zero_matrix(x_size, y_size);

    return matrix;
}

Matrix3D *
init_matrix3d(int x_size, int y_size, int z_size)
{
    Matrix3D *matrix = (Matrix3D *) malloc(sizeof(Matrix3D));

    matrix->z_size = z_size;
    matrix->y_size = y_size;
    matrix->x_size = x_size;
    matrix->data = init_zero_matrix3d(x_size, y_size, z_size);

    return matrix;
}

void
free_vector(Vector *vector)
{
    if (vector == NULL)
        return;
    
    free(vector->data);
    free(vector);
}

void
free_matrix(Matrix *matrix)
{
    if (matrix == NULL)
        return;

    for (int j = 0; j < matrix->y_size; j++)
        free(matrix->data[j]);
    
    free(matrix->data);
    free(matrix);
}

void
free_matrix3d(Matrix3D *matrix)
{
    if (matrix == NULL)
        return;

    for (int k = 0; k < matrix->z_size; k++)
    {
        for (int j = 0; j < matrix->y_size; j++)
            free(matrix->data[k][j]);
        free(matrix->data[k]);
    }
    
    free(matrix->data);
    free(matrix);
}

void
print_vector(Vector vector)
{
    for (int i = 0; i < vector.size; i++)
        printf("%.4f ", vector.data[i]);
    printf("\n");
}

void
print_matrix(Matrix matrix)
{
    for (int j = 0; j < matrix.y_size; j++)
    {
        for (int i = 0; i < matrix.x_size; i++)
            printf("%.4f ", matrix.data[j][i]);
        printf("\n");
    }
    printf("\n");
}

void
print_matrix3d(Matrix3D matrix)
{
    for (int k = 0; k < matrix.z_size; k++)
    {
        for (int j = 0; j < matrix.y_size; j++)
        {
            for (int i = 0; i < matrix.x_size; i++)
                printf("%.4f ", matrix.data[k][j][i]);
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

#ifdef USE_TIME_STRUCT
TimeVector *
init_time_vector(int time_size, int size)
{
    TimeVector *vector = (TimeVector *) malloc(sizeof(TimeVector));

    vector->time_size = time_size;
    vector->vector = malloc(sizeof(_vector) * time_size);
    for (int t = 0; t < time_size; t++)
        vector->vector[t] = init_vector(size);

    return vector;
}

TimeMatrix *
init_time_matrix(int time_size, int x_size, int y_size)
{
    TimeMatrix *matrix = (TimeMatrix *) malloc(sizeof(TimeMatrix));

    matrix->time_size = time_size;
    matrix->matrix = malloc(sizeof(_matrix) * time_size);
    for (int t = 0; t < time_size; t++)
        matrix->matrix[t] = init_matrix(x_size, y_size);

    return matrix;
}

TimeMatrix3D *
init_time_matrix3d(int time_size, int x_size, int y_size, int z_size)
{
    TimeMatrix3D *matrix = (TimeMatrix3D *) malloc(sizeof(TimeMatrix3D));
    
    matrix->time_size = time_size;
    matrix->matrix = malloc(sizeof(_matrix3D) * time_size);
    for (int t = 0; t < time_size; t++)
        matrix->matrix[t] = init_matrix3d(x_size, y_size, z_size);

    return matrix;
}

void
free_time_vector(TimeVector *vector)
{
    if (vector == NULL)
        return;

    for (int t = 0; t < vector->time_size; t++)
        free_vector(vector->vector[t]);
    free(vector->vector);
    free(vector);
}

void
free_time_matrix(TimeMatrix *matrix)
{
    if (matrix == NULL)
        return;

    for (int t = 0; t < matrix->time_size; t++)
        free_matrix(matrix->matrix[t]);
    free(matrix->matrix);
    free(matrix);
}

void
free_time_matrix3d(TimeMatrix3D *matrix)
{
    if (matrix == NULL)
        return;

    for (int t = 0; t < matrix->time_size; t++)
        free_matrix3d(matrix->matrix[t]);
    free(matrix->matrix);
    free(matrix);
}

void
print_time_vector(TimeVector vector, int step)
{
    if (step == 0)
    {
        printf("ERROR! Step for time cannot be equal zero!");
        return;
    }

    int vec_size = vector.size;
    for (int t = 0; t < vector.time_size; t += step)
    {
        for (int i = 0; i < vec_size; i++)
            printf("%.4f ", vector.vector[t][i]);
        printf("\n");
    }
    printf("\n");
}

void
print_time_matrix(TimeMatrix matrix, int step)
{
    if (step == 0)
    {
        printf("ERROR! Step for time cannot be equal zero!");
        return;
    }

    for (int t = 0; t < matrix.time_size; t += step)
    {
        for (int j = 0; j < matrix.y_size; j++)
        {
            for (int i = 0; i < matrix.x_size; i++)
                printf("%.4f ", matrix.matrix[t][j][i]);
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

void
print_time_matrix3d(TimeMatrix3D matrix, int step)
{
    if (step == 0)
    {
        printf("ERROR! Step for time cannot be equal zero!");
        return;
    }

    for (int t = 0; t < matrix.time_size; t += step)
    {
        for (int k = 0; k < matrix.z_size; k++)
        {
            for (int j = 0; j < matrix.y_size; j++)
            {
                for (int i = 0; i < matrix.x_size; i++)
                    printf("%.4f ", matrix.matrix[t][k][j][i]);
                printf("\n");
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}
#endif
