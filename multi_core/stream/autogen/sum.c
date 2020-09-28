#include <timing.h>
#include <likwid-marker.h>

#ifdef ACLE_VERSION
	#ifdef __ARM_FEATURE_SVE
		#include <arm_sve.h>
#GHOST_SUBST UNROLL_FACTOR #UNROLL_FACTOR_SUB#
	#else
		#error "SVE not supported by compiler"
	#endif /* __ARM_FEATURE_SVE */
#endif



double sum(
        double * restrict a,
        int N
        )
{
    double S, E;
    double sum = 0.0;

    S = getTimeStamp();


#ifdef ACLE_VERSION
    int pad = svcntd()*~UNROLL_FACTOR~;
    int N_round = ((int)(N/((double)pad)))*pad;
#endif

#pragma omp parallel
    {
        LIKWID_MARKER_START("SUM");

#ifdef ACLE_VERSION

	#GHOST_UNROLL#svfloat64_t sum_vec@ = svdup_f64(0);#UNROLL_FACTOR
	svbool_t pg = svptrue_b64();
#pragma omp for schedule(static) nowait
        for (int i=0; i<N_round; i+=svcntd()*~UNROLL_FACTOR~) {
		#GHOST_UNROLL#svfloat64_t ld_vec@ = svld1(pg, &(a[i+@*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#sum_vec@ = svadd_m(pg, sum_vec@, ld_vec@);#UNROLL_FACTOR
	}
	#GHOST_UNROLL#sum_vec0 = svadd_m(pg, sum_vec0, sum_vec~@+1~);#UNROLL_FACTOR-1
#pragma omp atomic
	sum += svaddv(svptrue_b64(), sum_vec0);

	//reminder loop
#pragma omp for reduction(+:sum) schedule(static)
        for (int i=N_round; i<N; i++) {
            sum += a[i];
        }
#else

#pragma omp for reduction(+:sum) schedule(static)
        for (int i=0; i<N; i++) {
            sum += a[i];
        }

#endif
        LIKWID_MARKER_STOP("SUM");
    }

    E = getTimeStamp();

    /* make the compiler think this makes actually sense */
    a[10] = sum;

    return E-S;
}
