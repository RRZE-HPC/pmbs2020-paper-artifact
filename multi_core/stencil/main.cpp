#include <iostream>
#include "timer.h"
#include "allocate.h"
#ifdef LIKWID_PERFMON
	#include "likwid.h"
#endif

@FWD_DECLARATION@

int main(int argc, char *argv[])
{
#ifdef LIKWID_PERFMON
	LIKWID_MARKER_INIT;
#endif
	double freq = 1.8*1e9;

	int N_start = 64;
	int N_end = 5000;

	if(argc <= 1)
	{
		printf("Usage : %s <unroll> <size(opt)>\n", argv[0]);
	}

	int unroll = atoi(argv[1]);
	if(argc > 2)
	{
		N_start = atoi(argv[2]);
		N_end = atoi(argv[2])+1;
	}
	printf("#%13s, %14s, %14s, %14s, %14s, %14s, %14s, %14s, %14s\n", "n_i", "n_j", "Size(kB)", "tot. time", "time/LUP(ms)", "cy/LUP", "cy/CL", "MLUPs", "rep");
//	int N = 1e8;
	for(int N_not_pad=N_start; N_not_pad<N_end; N_not_pad=N_not_pad*1.1)
	{
		double pad = unroll * 8;
		int N = ((int)(N_not_pad/pad))*pad;
		int n_j = N;
		int n_i = n_j; ///2.0;
		N = n_i*n_j;
		int N_alloc = N*1.2; //little extra alloc for avoiding extra predicates

#ifndef ALIGNED
		double *a = (double*) malloc(sizeof(double)*N_alloc);
		double *b = (double*) malloc(sizeof(double)*N_alloc);
#else
		double *a = (double *) allocate(1024, sizeof(double)*N_alloc);
		double *b = (double *) allocate(1024, sizeof(double)*N_alloc);
#endif

#pragma omp parallel for schedule(static)
		for(int i=0; i<n_i; ++i)
		{
			for(int j=0; j< n_j; ++j)
			{
				a[i*n_j+j] = 0.02*i;
				b[i*n_j+j] = 0.01*i;
			}
		}

		INIT_TIMER(load_warm);
		START_TIMER(load_warm);
#pragma omp parallel
		{
			for(int r=0; r<100; ++r)
			{
				@KERNEL@
			}
		}
		STOP_TIMER(load_warm);
		double warm_time = GET_TIMER(load_warm);
		int rep = 100*(0.5/warm_time);

		INIT_TIMER(load);
#ifdef LIKWID_PERFMON
		LIKWID_MARKER_START("load");
#endif
		START_TIMER(load);
#pragma omp parallel
		{
			for(int r=0; r<rep; ++r)
			{
				@KERNEL@
			}
		}
		STOP_TIMER(load);
#ifdef LIKWID_PERFMON
		LIKWID_MARKER_STOP("load");
#endif
		double time =  GET_TIMER(load);
		double numArrays = @N_ARRAYS@;
		double lup = (n_i-2)*(n_j-2);
		printf("%14d, %14d, %14.2f, %14.10f, %14.10f, %14.6f, %14.6f, %14f, %14d\n", n_i, n_j, numArrays*N*sizeof(double)/(1000.0), time, time*1e6/((double)N*rep), time*freq/((double)lup*rep), time*freq*8.0/((double)lup*rep), lup*rep*1e-6/(double)(time), rep);

		if(argc > 2)
		{
			printf("Rep = %d\n", rep);
			printf("n_i = %d\n", n_i);
			printf("n_j = %d\n", n_j);
			printf("lup = %f\n", lup);
			printf("N_size = %d\n", N);
			printf("Size = %f\n", numArrays*N*sizeof(double)/(1000.0));
			printf("Perf_cy = %f cy/CL\n", time*freq*8.0/((double)lup*rep));
			printf("MLUPS = %f \n", lup*rep*1e-6/(double)(time));
		}
		free(a);
		free(b);
	}

#ifdef LIKWID_PERFMON
	LIKWID_MARKER_CLOSE;
#endif
	return 0;
}
