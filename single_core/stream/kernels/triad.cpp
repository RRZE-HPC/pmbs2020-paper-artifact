#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif /* __ARM_FEATURE_SVE */

#GHOST_SUBST UNROLL_FACTOR #UNROLL_FACTOR_SUB#

void triad(double *a, double *b, double *c, double *d, int N)
{
	int64_t i = 0;
	int64_t n = N; //reduce N so I need only one pg
	svbool_t pg = svwhilelt_b64(i,n);
	double sum = 0;
	do
	{
		#GHOST_UNROLL#svfloat64_t b_vec@ = svld1(pg, &(b[i+~@~*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t c_vec@ = svld1(pg, &(c[i+~@~*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t d_vec@ = svld1(pg, &(d[i+~@~*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t a_vec@ = svmad_m(pg, d_vec@, c_vec@, b_vec@);#UNROLL_FACTOR 
		#GHOST_UNROLL#svst1(pg, &(a[i+~@~*svcntd()]), a_vec@);#UNROLL_FACTOR

		i += svcntd()*~UNROLL_FACTOR~;
		pg = svwhilelt_b64(i, n);
	} while(svptest_any(svptrue_b64(), pg));

	return;
}
