#ifndef ARRAY_UTILS
#define ARRAY_UTILS
float *init_zero_1d(int size);

float **init_zero_2d(int x_size, int y_size);

float ***init_zero_3d(int x_size, int y_size, int z_size);

void free_1d_array(float *array);

void free_2d_array(float **array, int M);

void free_3d_array(float ***array, int K, int M);

void print_array_1d(float *array, int size);

void print_array_2d(float **array, int M, int N);

void print_array_3d(float ***array, int K, int M, int N);

#endif