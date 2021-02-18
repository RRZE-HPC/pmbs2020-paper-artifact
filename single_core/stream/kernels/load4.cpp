#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif /* __ARM_FEATURE_SVE */

#GHOST_SUBST UNROLL_FACTOR #UNROLL_FACTOR_SUB#

double load4(double *a, double *b,double *c, double *d, int N)
{
	int64_t i = 0;
	int64_t n = N; //reduce N so I need only one pg
	svbool_t pg = svwhilelt_b64(i,n);
	double sum = 0;
	#GHOST_UNROLL#float64_t sum@ = 0;#UNROLL_FACTOR
	do
	{
		#GHOST_UNROLL#svfloat64_t ld1_vec@ = svld1(pg, &(a[i+@*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t ld2_vec@ = svld1(pg, &(b[i+@*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t ld3_vec@ = svld1(pg, &(c[i+@*svcntd()]));#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t ld4_vec@ = svld1(pg, &(d[i+@*svcntd()]));#UNROLL_FACTOR
	
		#GHOST_UNROLL#svfloat64_t sum1_vec@ = svadd_m(pg, ld1_vec@, ld2_vec@);#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t sum2_vec@ = svadd_m(pg, ld3_vec@, ld4_vec@);#UNROLL_FACTOR
		#GHOST_UNROLL#svfloat64_t sum_vec@ = svadd_m(pg, sum1_vec@, sum2_vec@);#UNROLL_FACTOR
		#GHOST_UNROLL#sum@ += svaddv(svptrue_b64(), sum_vec@);#UNROLL_FACTOR

		i += svcntd()*~UNROLL_FACTOR~;
		pg = svwhilelt_b64(i, n);
	} while(svptest_any(svptrue_b64(), pg));

	#GHOST_UNROLL#sum += sum@;#UNROLL_FACTOR
	return sum;
}
