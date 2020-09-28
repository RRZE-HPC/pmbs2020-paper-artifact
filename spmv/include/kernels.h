#ifndef _KERNELS_H
#define _KERNELS_H

#include "sparsemat.h"
#include "densemat.h"
#include <stdio.h>
#include "macros.h"

void spmv(densemat* b, sparsemat* mat, densemat* x);
void spmv_a64fx(densemat* b, sparsemat* mat, densemat* x);
void spmv_a64fx_sell32(densemat* b, sparsemat* mat, densemat* x);

#endif
