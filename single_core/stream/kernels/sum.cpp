#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif /* __ARM_FEATURE_SVE */

#GHOST_SUBST UNROLL_FACTOR #UNROLL_FACTOR_SUB#

double sum(double *a, int N)
{
	int64_t i = 0;
	int64_t n = N; //reduce N so I need only one pg
	svbool_t pg = svwhilelt_b64(i,n);
	double sum = 0;
	#GHOST_UNROLL#svfloat64_t sum_vec@ = svdup_f64_x(svptrue_b64(), 0);#UNROLL_FACTOR
	do
	{
		#GHOST_UNROLL#svfloat64_t ld_vec@ = svld1(pg, &(a[i+@*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#sum_vec@ = svadd_m(pg, sum_vec@, ld_vec@);#UNROLL_FACTOR

		i += svcntd()*~UNROLL_FACTOR~;
		pg = svwhilelt_b64(i, n);
	} while(svptest_any(svptrue_b64(), pg));

	#GHOST_UNROLL#float64_t sum@ = svaddv(svptrue_b64(), sum_vec@);#UNROLL_FACTOR
	#GHOST_UNROLL#sum += sum@;#UNROLL_FACTOR
	return sum;
}
