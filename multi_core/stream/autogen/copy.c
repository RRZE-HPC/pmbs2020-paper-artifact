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


double copy(
        double * restrict a,
        double * restrict b,
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
        LIKWID_MARKER_START("COPY");
#ifdef ACLE_VERSION

	svbool_t pg = svptrue_b64();
#pragma omp for schedule(static)
        for (int i=0; i<N_round; i+=svcntd()*~UNROLL_FACTOR~) {
		#GHOST_UNROLL#svfloat64_t b_vec@ = svld1(pg, &(b[i+@*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#svst1(pg, &(a[i+~@~*svcntd()]), b_vec@);#UNROLL_FACTOR
	}

#pragma omp for schedule(static)
        for (int i=N_round; i<N; i++) {
            a[i] = b[i];
        }
#else


#pragma omp for
        for (int i=0; i<N; i++) {
            a[i] = b[i];
        }
#endif
        LIKWID_MARKER_STOP("COPY");
    }
    E = getTimeStamp();

    return E-S;
}
