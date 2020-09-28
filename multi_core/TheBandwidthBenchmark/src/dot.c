/*
 * =======================================================================================
 *
 *      Author:   Jan Eitzinger (je), jan.eitzinger@fau.de
 *      Copyright (c) 2019 RRZE, University Erlangen-Nuremberg
 *
 *      Permission is hereby granted, free of charge, to any person obtaining a copy
 *      of this software and associated documentation files (the "Software"), to deal
 *      in the Software without restriction, including without limitation the rights
 *      to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *      copies of the Software, and to permit persons to whom the Software is
 *      furnished to do so, subject to the following conditions:
 *
 *      The above copyright notice and this permission notice shall be included in all
 *      copies or substantial portions of the Software.
 *
 *      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *      IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *      FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *      AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *      LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *      OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *      SOFTWARE.
 *
 * =======================================================================================
 */

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
    int pad = svcntd()*4;
    int N_round = ((int)(N/((double)pad)))*pad;
#endif

#pragma omp parallel
    {
        LIKWID_MARKER_START("DOT");

#ifdef ACLE_VERSION

	svfloat64_t sum_vec0 = svdup_f64(0);
	svfloat64_t sum_vec1 = svdup_f64(0);
	svfloat64_t sum_vec2 = svdup_f64(0);
	svfloat64_t sum_vec3 = svdup_f64(0);
	svbool_t pg = svptrue_b64();
#pragma omp for schedule(static) nowait
        for (int i=0; i<N_round; i+=svcntd()*4) {
		svfloat64_t a_vec0 = svld1(pg, &(a[i+0*svcntd()]));
		svfloat64_t a_vec1 = svld1(pg, &(a[i+1*svcntd()]));
		svfloat64_t a_vec2 = svld1(pg, &(a[i+2*svcntd()]));
		svfloat64_t a_vec3 = svld1(pg, &(a[i+3*svcntd()]));
		svfloat64_t b_vec0 = svld1(pg, &(b[i+0*svcntd()]));
		svfloat64_t b_vec1 = svld1(pg, &(b[i+1*svcntd()]));
		svfloat64_t b_vec2 = svld1(pg, &(b[i+2*svcntd()]));
		svfloat64_t b_vec3 = svld1(pg, &(b[i+3*svcntd()]));
		sum_vec0 = svmla_m(pg, sum_vec0, a_vec0, b_vec0);
		sum_vec1 = svmla_m(pg, sum_vec1, a_vec1, b_vec1);
		sum_vec2 = svmla_m(pg, sum_vec2, a_vec2, b_vec2);
		sum_vec3 = svmla_m(pg, sum_vec3, a_vec3, b_vec3);
	}
	sum_vec0 = svadd_m(pg, sum_vec0, sum_vec1);
	sum_vec0 = svadd_m(pg, sum_vec0, sum_vec2);
	sum_vec0 = svadd_m(pg, sum_vec0, sum_vec3);
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
