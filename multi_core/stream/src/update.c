#include <timing.h>
#include <likwid-marker.h>

#ifdef ACLE_VERSION
	#ifdef __ARM_FEATURE_SVE
		#include <arm_sve.h>
	#else
		#error "SVE not supported by compiler"
	#endif /* __ARM_FEATURE_SVE */
#endif


double update(
        double * restrict a,
        double scalar,
        int N
        )
{
    double S, E;

    S = getTimeStamp();

#ifdef ACLE_VERSION
    int pad = svcntd()*1;
    int N_round = ((int)(N/((double)pad)))*pad;
#endif

#pragma omp parallel
    {
        LIKWID_MARKER_START("UPDATE");

#ifdef ACLE_VERSION

	svbool_t pg = svptrue_b64();
	svfloat64_t scalar_vec  = svdup_f64_x(pg, scalar);
#pragma omp for schedule(static)
        for (int i=0; i<N_round; i+=svcntd()*1) {
		svfloat64_t a_vec0 = svld1(pg, &(a[i+0*svcntd()]));
		a_vec0 = svmul_m(pg, scalar_vec, a_vec0);
		svst1(pg, &(a[i+0*svcntd()]), a_vec0);
	}

#pragma omp for schedule(static)
        for (int i=N_round; i<N; i++) {
            a[i] = a[i] * scalar;
        }

#else


#pragma omp for schedule(static)
        for (int i=0; i<N; i++) {
            a[i] = a[i] * scalar;
        }
#endif
        LIKWID_MARKER_STOP("UPDATE");
    }
    E = getTimeStamp();

    return E-S;
}
