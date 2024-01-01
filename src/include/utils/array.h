#ifndef UTILS_H
#define UTILS_H

typedef float* Vector;
typedef float** Matrix; 
typedef float*** Matrix3D; 


Vector init_zero_vector(int size);
Matrix init_zero_matrix(int x_size, int y_size);
Matrix3D init_zero_matrix3d(int x_size, int y_size, int z_size);

void free_vector(Vector vector);
void free_matrix(Matrix matrix);
void free_matrix3d(Matrix3D array);

void print_vector(Vector vector, int size);
void print_matrix(Matrix matrix);
void print_matrix3d(Matrix3D array);

#endif /* UTILS_H */