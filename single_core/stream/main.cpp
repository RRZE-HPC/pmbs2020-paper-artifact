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
	double freq = 2.2*1e9;

	int N_start = 200;
	int N_end = 56000000;

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
	printf("#%13s, %14s, %14s, %14s, %14s, %14s, %14s\n", "N", "Size(kB)", "tot. time", "time/LUP(ms)", "cy/LUP", "cy/CL", "rep");
//	int N = 1e8;
	for(int N_not_pad=N_start; N_not_pad<N_end; N_not_pad=N_not_pad*1.2)
	{
		double pad = unroll * 8;
		int N = ((int)(N_not_pad/pad))*pad;
		int N_alloc = N*1.2; //little extra alloc for avoiding extra predicates

#ifndef ALIGNED
		double *a = (double*) malloc(sizeof(double)*N_alloc);
		double *b = (double*) malloc(sizeof(double)*N_alloc);
		double *c = (double*) malloc(sizeof(double)*N_alloc);
		double *d = (double*) malloc(sizeof(double)*N_alloc);
#elif defined(THPALLOC)
		double *a = (double *) hp_allocate(sizeof(double)*N_alloc);
		double *b = (double *) hp_allocate(sizeof(double)*N_alloc);
		double *c = (double *) hp_allocate(sizeof(double)*N_alloc);
		double *d = (double *) hp_allocate(sizeof(double)*N_alloc);
#else
		double *a = (double *) allocate(1024, sizeof(double)*N_alloc);
		double *b = (double *) allocate(1024, sizeof(double)*N_alloc);
		double *c = (double *) allocate(1024, sizeof(double)*N_alloc);
		double *d = (double *) allocate(1024, sizeof(double)*N_alloc);
#endif
		for(int i=0; i<N_alloc; ++i)
		{
			a[i] = 0.04*i;
			b[i] = 0.03*i;
			c[i] = 0.02*i;
			d[i] = 0.01*i;
		}

		INIT_TIMER(load_warm);
		START_TIMER(load_warm);
		for(int r=0; r<100; ++r)
		{
			@KERNEL@
		}
		STOP_TIMER(load_warm);
		double warm_time = GET_TIMER(load_warm);
		int rep = 100*(0.5/warm_time);

		INIT_TIMER(load);
#ifdef LIKWID_PERFMON
		LIKWID_MARKER_START("load");
#endif
		START_TIMER(load);
		for(int r=0; r<rep; ++r)
		{
			@KERNEL@
		}
		STOP_TIMER(load);
#ifdef LIKWID_PERFMON
		LIKWID_MARKER_STOP("load");
#endif
		double time =  GET_TIMER(load);
		double numArrays = @N_ARRAYS@;
		printf("%14d, %14.2f, %14.10f, %14.10f, %14.6f, %14.6f, %14d\n", N, numArrays*N*sizeof(double)/(1000.0), time, time*1e6/((double)N*rep), time*freq/((double)N*rep), time*freq*8.0/((double)N*rep), rep);

		if(argc > 2)
		{
			printf("Rep = %d\n", rep);
			printf("N_size = %d\n", N);
			printf("Size = %f\n", numArrays*N*sizeof(double)/(1000.0));
			printf("Perf_cy = %f cy/CL\n", time*freq*8.0/((double)N*rep));
		}
#ifndef THPALLOC
		free(a);
		free(b);
		free(c);
		free(d);
#else
        hp_free(a, sizeof(double)*N_alloc);
        hp_free(b, sizeof(double)*N_alloc);
        hp_free(c, sizeof(double)*N_alloc);
        hp_free(d, sizeof(double)*N_alloc);
#endif
	}

#ifdef LIKWID_PERFMON
	LIKWID_MARKER_CLOSE;
#endif
	return 0;
}
