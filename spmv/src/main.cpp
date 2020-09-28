#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include "mmio.h"
#include "time.h"
#include <cctype>
#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif
#include "parse.h"
#include "sparsemat.h"
#include "densemat.h"
#include "kernels.h"
#include "timer.h"
#include "spMV_NLPE.h"
#include <string.h>

//#include <mpi.h>
void capitalize(char* beg)
{
    int i = 0;
    while(beg[i])
    {
        beg[i] = toupper(beg[i]);
        ++i;
    }
}

//overload new and delete for alignement
void * operator new(size_t bytesize)
{
	//printf("Overloading new operator with size: %lu\n", bytesize);
	int errorCode;
	void* ptr;
	int alignment = 1024;
	errorCode =  posix_memalign(&ptr, alignment, bytesize);

	if (errorCode) {
		if (errorCode == EINVAL) {
			fprintf(stderr,
					"Error: Alignment parameter is not a power of two\n");
			exit(EXIT_FAILURE);
		}
		if (errorCode == ENOMEM) {
			fprintf(stderr,
					"Error: Insufficient memory to fulfill the request\n");
			exit(EXIT_FAILURE);
		}
	}

	if (ptr == NULL) {
		fprintf(stderr, "Error: posix_memalign failed!\n");
		exit(EXIT_FAILURE);
	}

	return ptr;
}

void operator delete(void * p)
{
	//printf("Overloading delete operator\n");
	free(p);
}

#define EXCHANGE(_lhs_, _rhs_)\
{\
	double* tmp = _rhs_->val;\
	_rhs_->val = _lhs_->val;\
	_lhs_->val = tmp;\
}\



#ifdef LIKWID_PERFMON

#define PERF_RUN(kernel, flopPerNnz, _lhs_)\
{\
    /*  int myRank = 0;\
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);*/\
    double time = 0;\
    double nnz_update = ((double)mat->nnz)*iter*1e-9;\
    sleep(1);\
    _Pragma("omp parallel")\
    {\
        LIKWID_MARKER_START(#kernel);\
    }\
    INIT_TIMER(kernel);\
    START_TIMER(kernel);\
    /*_Pragma("omp parallel")*/\
    {\
        for(int it=0; it<iter; ++it)\
        {\
            kernel(_lhs_, mat, x);\
            EXCHANGE(_lhs_, x);\
        }\
    }\
    STOP_TIMER(kernel);\
    time = GET_TIMER(kernel);\
    _Pragma("omp parallel")\
    {\
        LIKWID_MARKER_STOP(#kernel);\
    }\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
    /*   printf("%10s(%d) : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, myRank, flopPerNnz*nnz_update/(time), time);*/\
    printf("%10s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);\
    free(capsKernel);\
}\

#else

#define PERF_RUN(kernel, flopPerNnz, _lhs_)\
{\
    /*  int myRank = 0;\
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);*/\
    double time = 0;\
    double nnz_update = ((double)mat->nnz)*iter*1e-9;\
    sleep(1);\
    INIT_TIMER(kernel);\
    START_TIMER(kernel);\
    /*_Pragma("omp parallel")*/\
    {\
        for(int it=0; it<iter; ++it)\
        {\
            kernel(_lhs_, mat, x);\
            EXCHANGE(_lhs_, x);\
        }\
    }\
    STOP_TIMER(kernel);\
    time = GET_TIMER(kernel);\
    char* capsKernel;\
    asprintf(&capsKernel, "%s", #kernel);\
    capitalize(capsKernel);\
    /*   printf("%10s(%d) : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, myRank, flopPerNnz*nnz_update/(time), time);*/\
    printf("%10s : %8.4f GFlop/s ; Time = %8.5f s\n", capsKernel, flopPerNnz*nnz_update/(time), time);\
    free(capsKernel);\
}\

#endif

int main(int argc, char * argv[])
{

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_INIT;
#endif
    int err;
    parser param;
    if(!param.parse_arg(argc, argv))
    {
        printf("Error in reading parameters\n");
    }
    sparsemat* mat = new sparsemat;


    if(param.mat_file)
    {
	    printf("Reading matrix file\n");
	    if(!mat->readFile(param.mat_file))
	    {
		    printf("Error in reading sparse matrix file\n");
	    }
    }
    else
    {
        int S = param.hpcg_size;
        printf("Creating HPCG with dimension %d\n", S);
        mat->generateHPCG(S, S, S);
    }
    if(param.RCM_flag)
    {
        mat->doRCMPermute();
    }

    if(param.chunkHeight>1)
    {
        printf("Constructing SELL-C-sigma with C = %d, sigma = %d\n", param.chunkHeight, param.sigma);
        mat->constructSellCSigma(param.chunkHeight, param.sigma);
    }

    int NROWS = mat->nrows;
    if(param.chunkHeight>1)
    { 
        NROWS = mat->nchunks*mat->C;
    }
    if(param.chunkHeight == 1)
    {
        printf("nrows_Sell_C = %d, nnz_Sell_C = %d, nnzr_Sell_C = %f\n", NROWS, mat->nnz, mat->nnz/(double)NROWS);
        printf("nrows_CRS = %d, nnz_CRS = %d, nnzr_CRS = %f\n", NROWS, mat->nnz, mat->nnz/(double)mat->nrows);
        printf("eff = %f\n", 1.0);
    }
    else
    {
        printf("nrows_Sell_C = %d, nnz_Sell_C = %d, nnzr_Sell_C = %f\n", NROWS, mat->nnzSellC, mat->nnzSellC/(double)NROWS);
        printf("nrows_CRS = %d, nnz_CRS = %d, nnzr_CRS = %f\n", mat->nrows, mat->nnz, mat->nnz/(double)mat->nrows);
        printf("eff = %f\n", mat->nnz/(double)(mat->nnzSellC));
    }

    densemat *x = NULL, *b_a64fx;
    x=new densemat(NROWS, 1);
    b_a64fx=new densemat(NROWS, 1);

    double xinit = 1;
    x->setVal(xinit);
    //x->setVal(1);
    b_a64fx->setRand();

    //find iter
    int iter = 10;
    INIT_TIMER(warm);
    START_TIMER(warm);
    //#pragma omp parallel
    {
        for(int it=0; it<iter; ++it)
        {
            spmv_a64fx(b_a64fx, mat, x);
            EXCHANGE(b_a64fx, x);
        }
    }

    x->setVal(xinit);
    STOP_TIMER(warm);
    double warm_time = GET_TIMER(warm);
    iter = 1*iter/warm_time;
    iter = (iter%2 == 0)?iter:iter+1; //make it even, so validation works
    printf("Iter = %d\n", iter);
    //This macro times and reports performance
    PERF_RUN(spmv_a64fx,2, b_a64fx);

    if(param.validate)
    {
        densemat *b=new densemat(NROWS, 1);
        b->setRand();
        x->setVal(xinit);
        PERF_RUN(spmv,2, b);
        checkEqual(b,b_a64fx, xinit);

	delete b;
    }

    delete mat;
    delete b_a64fx;
    delete x;

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif
}
