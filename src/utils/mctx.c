#include <stdlib.h>

void *
conv_malloc(int size)
{
    void *res = malloc(sizeof(size));

    return res;
}