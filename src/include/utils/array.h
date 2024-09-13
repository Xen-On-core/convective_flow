#ifndef ARRAYS_H
#define ARRAYS_H

typedef double* _vector;
typedef double** _matrix;
typedef double*** _matrix3D;

typedef struct Vector
{
    _vector data;
    int size;
} Vector;

typedef struct Matrix
{
    _matrix data;
    int x_size;
    int y_size;
} Matrix;

typedef struct Matrix3D
{
    _matrix3D data;
    int x_size;
    int y_size;
    int z_size;
} Matrix3D;

Vector *init_vector(int size);
Matrix *init_matrix(int x_size, int y_size);
Matrix3D *init_matrix3d(int x_size, int y_size, int z_size);

void free_vector(Vector *vector);
void free_matrix(Matrix *matrix);
void free_matrix3d(Matrix3D *matrix);

void print_vector(Vector vector);
void print_matrix(Matrix matrix);
void print_matrix3d(Matrix3D matrix);

    #ifdef USE_TIME_STRUCT
typedef struct TimeVector
{
    Vector *vector;
    int size;
    int time_size;
} TimeVector;

typedef struct TimeMatrix
{
    _matrix *matrix;
    int x_size;
    int y_size;
    int time_size;
} TimeMatrix;

typedef struct TimeMatrix3D
{
    _matrix3D *matrix;
    int x_size;
    int y_size;
    int z_size;
    int time_size;
} TimeMatrix3D;

TimeVector *init_time_vector(int size, int time_size);
TimeMatrix *init_time_matrix(int x_size, int y_size, int time_size);
TimeMatrix3D *init_time_matrix3d(int x_size, int y_size, int z_size, int time_size);

void free_time_vector(TimeVector *vector);
void free_time_matrix(TimeMatrix *matrix);
void free_time_matrix3d(TimeMatrix3D *matrix);

void print_time_vector(TimeVector vector, int step);
void print_time_matrix(TimeMatrix matrix, int step);
void print_time_matrix3d(TimeMatrix3D matrix, int step);
    #endif

#endif /* ARRAYS_H */