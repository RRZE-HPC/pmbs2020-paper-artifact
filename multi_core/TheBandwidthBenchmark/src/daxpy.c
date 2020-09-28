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
    int pad = svcntd()*4;
    int N_round = ((int)(N/((double)pad)))*pad;
#endif

#pragma omp parallel
    {
        LIKWID_MARKER_START("DAXPY");
#ifdef ACLE_VERSION

	svbool_t pg = svptrue_b64();
	svfloat64_t scalar_vec  = svdup_f64_x(pg, scalar);
#pragma omp for schedule(static)
        for (int i=0; i<N_round; i+=svcntd()*4) {
		svfloat64_t a_vec0 = svld1(pg, &(a[i+0*svcntd()]));
		svfloat64_t a_vec1 = svld1(pg, &(a[i+1*svcntd()]));
		svfloat64_t a_vec2 = svld1(pg, &(a[i+2*svcntd()]));
		svfloat64_t a_vec3 = svld1(pg, &(a[i+3*svcntd()]));
		svfloat64_t b_vec0 = svld1(pg, &(b[i+0*svcntd()]));
		svfloat64_t b_vec1 = svld1(pg, &(b[i+1*svcntd()]));
		svfloat64_t b_vec2 = svld1(pg, &(b[i+2*svcntd()]));
		svfloat64_t b_vec3 = svld1(pg, &(b[i+3*svcntd()]));
		a_vec0 = svmad_m(pg, scalar_vec, b_vec0, a_vec0);
		a_vec1 = svmad_m(pg, scalar_vec, b_vec1, a_vec1);
		a_vec2 = svmad_m(pg, scalar_vec, b_vec2, a_vec2);
		a_vec3 = svmad_m(pg, scalar_vec, b_vec3, a_vec3);
		svst1(pg, &(a[i+0*svcntd()]), a_vec0);
		svst1(pg, &(a[i+1*svcntd()]), a_vec1);
		svst1(pg, &(a[i+2*svcntd()]), a_vec2);
		svst1(pg, &(a[i+3*svcntd()]), a_vec3);
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
