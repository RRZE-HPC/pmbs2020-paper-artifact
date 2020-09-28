#include "densemat.h"
#include "math.h"
#include <cstdio>
#include <omp.h>

//ncols is for nodes
densemat::densemat(int nrows_, int ncols_):nrows(nrows_), ncols(ncols_)
{
	val = new double[nrows*ncols];

	if(ncols != 1)
	{
		printf("Error: More than 1 columnn not yet supported\n");
	}
#pragma omp parallel for schedule(static)
	for(int i=0; i<nrows; ++i)
	{
		val[i] = 0.0;
	}
}

densemat::~densemat()
{
    if(val)
    {
        delete[] val;
    }
}

void densemat::setVal(double value)
{
    for(int j=0; j<ncols; ++j)
    {
#pragma omp parallel for schedule(static)
        for(int i=0; i<nrows; ++i)
        {
            val[j*nrows+i] = value;
        }
    }
}
/*
   void densemat::setFn(std::function<double(int)> fn)
   {
#pragma omp parallel for
for(int i=0; i<nrows; ++i)
{
val[i] = fn(i);
}
}*/

void densemat::setFn(std::function<double(void)> fn)
{
    for(int j=0; j<ncols; ++j)
    {
#pragma omp parallel for schedule(static)
        for(int i=0; i<nrows; ++i)
        {
            val[j*nrows+i] = fn();
        }
    }
}

void densemat::setRand()
{
    setFn(rand);
}

bool checkEqual(const densemat* lhs, const densemat* rhs, double tol)
{
    if((lhs->nrows != rhs->nrows) && (lhs->ncols != rhs->ncols))
    {
        printf("Densemat dimension differs\n");
        return false;
    }

    int nrows = lhs->nrows;
    int ncols = lhs->ncols;

    for(int col=0; col<ncols; ++col)
    {
        for(int row=0; row<nrows; ++row)
        {
            if( fabs((lhs->val[col*nrows+row]-rhs->val[col*nrows+row])/lhs->val[col*nrows+row]) > tol )
            {
                printf("Densemat deviation @ idx %d lhs = %f, rhs = %f\n", row, lhs->val[col*nrows+row], rhs->val[col*nrows+row]);
                return false;
            }
        }
    }

    return true;
}
