#ifndef ARRAYS_H
#define ARRAYS_H

typedef double* Vector;
typedef double** Matrix;
typedef double*** Matrix3D;

Vector init_zero_vector(int size);
Matrix init_zero_matrix(int x_size, int y_size);
Matrix3D init_zero_matrix3d(int x_size, int y_size, int z_size);

void free_vector(Vector vector);
void free_matrix(Matrix matrix);
void free_matrix3d(Matrix3D matrix);

void print_vector(Vector vector);
void print_matrix(Matrix matrix);
void print_matrix3d(Matrix3D matrix);

typedef Vector* TimeVector;
typedef Matrix* TimeMatrix;
typedef Matrix3D* TimeMatrix3D;

TimeVector init_zero_time_vector(int size, int time_size);
TimeMatrix init_zero_time_matrix(int x_size, int y_size, int time_size);
TimeMatrix3D init_zero_time_matrix3d(int x_size, int y_size, int z_size, int time_size);

void free_time_vector(TimeVector vector);
void free_time_matrix(TimeMatrix matrix);
void free_time_matrix3d(TimeMatrix3D matrix);

void print_time_vector(TimeVector vector, int step);
void print_time_matrix(TimeMatrix matrix, int step);
void print_time_matrix3d(TimeMatrix3D matrix, int step);

#endif /* ARRAYS_H */