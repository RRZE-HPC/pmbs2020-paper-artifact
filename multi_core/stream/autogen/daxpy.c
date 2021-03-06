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


double daxpy(
        double * restrict a,
        double * restrict b,
        double scalar,
        int N
        )
{
    double S, E;

    S = getTimeStamp();

#ifdef ACLE_VERSION
    int pad = svcntd()*~UNROLL_FACTOR~;
    int N_round = ((int)(N/((double)pad)))*pad;
#endif

#pragma omp parallel
    {
        LIKWID_MARKER_START("DAXPY");
#ifdef ACLE_VERSION

	svbool_t pg = svptrue_b64();
	svfloat64_t scalar_vec  = svdup_f64_x(pg, scalar);
#pragma omp for schedule(static)
        for (int i=0; i<N_round; i+=svcntd()*~UNROLL_FACTOR~) {
		#GHOST_UNROLL#svfloat64_t a_vec@ = svld1(pg, &(a[i+@*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t b_vec@ = svld1(pg, &(b[i+@*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#a_vec@ = svmad_m(pg, scalar_vec, b_vec@, a_vec@);#UNROLL_FACTOR
		#GHOST_UNROLL#svst1(pg, &(a[i+~@~*svcntd()]), a_vec@);#UNROLL_FACTOR
	}

#pragma omp for schedule(static)
        for (int i=N_round; i<N; i++) {
            a[i] = a[i] + scalar * b[i];
        }

#else


#pragma omp for schedule(static)
        for (int i=0; i<N; i++) {
            a[i] = a[i] + scalar * b[i];
        }
#endif
	LIKWID_MARKER_STOP("DAXPY");
    }
    E = getTimeStamp();

    return E-S;
}
