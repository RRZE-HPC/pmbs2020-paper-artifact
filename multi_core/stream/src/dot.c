#include <timing.h>
#include <likwid-marker.h>

#ifdef ACLE_VERSION
	#ifdef __ARM_FEATURE_SVE
		#include <arm_sve.h>
	#else
		#error "SVE not supported by compiler"
	#endif /* __ARM_FEATURE_SVE */
#endif



double dot(
        double * restrict a,
        double * restrict b,
        int N
        )
{
    double S, E;
    double sum = 0.0;

    S = getTimeStamp();


#ifdef ACLE_VERSION
    int pad = svcntd()*1;
    int N_round = ((int)(N/((double)pad)))*pad;
#endif

#pragma omp parallel
    {
        LIKWID_MARKER_START("DOT");

#ifdef ACLE_VERSION

	svfloat64_t sum_vec0 = svdup_f64(0);
	svbool_t pg = svptrue_b64();
#pragma omp for schedule(static) nowait
        for (int i=0; i<N_round; i+=svcntd()*1) {
		svfloat64_t a_vec0 = svld1(pg, &(a[i+0*svcntd()]));
		svfloat64_t b_vec0 = svld1(pg, &(b[i+0*svcntd()]));
		sum_vec0 = svmla_m(pg, sum_vec0, a_vec0, b_vec0);
	}
#pragma omp atomic
	sum += svaddv(svptrue_b64(), sum_vec0);

	//reminder loop
#pragma omp for reduction(+:sum) schedule(static)
        for (int i=N_round; i<N; i++) {
            sum += a[i]*b[i];
        }
#else

#pragma omp for reduction(+:sum) schedule(static)
        for (int i=0; i<N; i++) {
            sum += a[i]*b[i];
        }

#endif
        LIKWID_MARKER_STOP("DOT");
    }

    E = getTimeStamp();

    /* make the compiler think this makes actually sense */
    a[10] = sum;

    return E-S;
}
