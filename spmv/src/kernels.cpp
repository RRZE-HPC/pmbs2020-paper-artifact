#include "kernels.h"
#include "print.h"
#include <omp.h>

inline void spmv_csr_general(densemat* b, sparsemat* mat, densemat* x)
{

    INFO_PRINT("CRS General");
#pragma omp parallel for schedule(static)
    for(int row=0; row<mat->nrows; ++row)
    {
        double tmp = 0;
        for(int idx=mat->rowPtr[row]; idx<mat->rowPtr[row+1]; ++idx)
        {
            tmp += mat->val[idx]*x->val[mat->col[idx]];
        }
        b->val[row] = tmp;
    }
}

inline void spmv_sellC_general(densemat* b, sparsemat* mat, densemat* x)
{
    INFO_PRINT("SELL-C General");
    const int C = mat->C;
    const int  P = mat->P;

#pragma omp parallel for schedule(static)
    for(int chunk=0; chunk<mat->nchunks; ++chunk)
    {
        for(int rowInChunk=0; rowInChunk<C; ++rowInChunk)
        {
            b->val[chunk*C+rowInChunk] = 0;
        }
        for(int j=0; j<mat->chunkLen[chunk]; j=j+P)
        {
            int idx = mat->chunkPtr[chunk]+j*C;
            for(int rowInChunk=0; rowInChunk<C; ++rowInChunk)
            {
                b->val[chunk*C+rowInChunk] += mat->valSellC[idx+rowInChunk]*x->val[mat->colSellC[idx+rowInChunk]];
            }
        }
    }
}


//b=A*x
void spmv(densemat* b, sparsemat* mat, densemat* x)
{
	if(mat->C==1)
	{
		spmv_csr_general(b, mat, x);
	}
	else
	{
		spmv_sellC_general(b, mat, x);
	}
}

void spmv_a64fx(densemat* b, sparsemat* mat, densemat* x)
{
    if(mat->C==1)
    {
	    spmv(b, mat, x);
    }
    else
    {
	    switch(mat->C)
	    {
		    case 32:
			    spmv_a64fx_sell32(b, mat, x);
			    break;
		    default:
			    WARNING_PRINT("Requested SELL-%d kernel with lb not found, calling fallback\n", mat->C);
			    spmv(b, mat, x);
			    break;
	    }
    }

}


