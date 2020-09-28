#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif /* __ARM_FEATURE_SVE */

#GHOST_SUBST UNROLL_FACTOR #UNROLL_FACTOR_SUB#

void stencil_2d5pt(double *a, double *b, int n_i, int n_j)
{

	for(int i=1; i<(n_i-1); ++i)
	{
		int64_t j = 0;
		int64_t n_j_64 = n_j;
		svbool_t pg = svwhilelt_b64(j,n_j_64);


		do
		{
			#GHOST_UNROLL#svfloat64_t bottom@ = svld1(pg, &(b[(i-1)*n_j+j+~@~*svcntd()]));#UNROLL_FACTOR
			#GHOST_UNROLL#svfloat64_t left@ = svld1(pg, &(b[i*n_j+(j-1)+~@~*svcntd()]));#UNROLL_FACTOR
			#GHOST_UNROLL#svfloat64_t center@ = svld1(pg, &(b[i*n_j+(j)+~@~*svcntd()]));#UNROLL_FACTOR
			#GHOST_UNROLL#svfloat64_t right@ = svld1(pg, &(b[i*n_j+(j+1)+~@~*svcntd()]));#UNROLL_FACTOR
			#GHOST_UNROLL#svfloat64_t top@ = svld1(pg, &(b[(i+1)*n_j+j+~@~*svcntd()]));#UNROLL_FACTOR

			#GHOST_UNROLL#bottom@ = svadd_m(pg, bottom@, left@);#UNROLL_FACTOR
			#GHOST_UNROLL#center@ = svadd_m(pg, center@, right@);#UNROLL_FACTOR
			#GHOST_UNROLL#top@ = svadd_m(pg, bottom@, top@);#UNROLL_FACTOR
			#GHOST_UNROLL#center@ = svadd_m(pg, center@, top@);#UNROLL_FACTOR

			#GHOST_UNROLL#svst1(pg, &(a[(i)*n_j+j+~@~*svcntd()]), center@);#UNROLL_FACTOR

			j += svcntd()*~UNROLL_FACTOR~;
			pg = svwhilelt_b64(j, n_j_64);
		} while(svptest_any(svptrue_b64(), pg));
	}
	return;
}
