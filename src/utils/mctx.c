#include <stdlib.h>
#include "utils/mctx.h"

void *
cnvalloc(int size)
{
    void *res = malloc(size);
    if (res == NULL)
        alloc_error = 1;

    return res;
}
