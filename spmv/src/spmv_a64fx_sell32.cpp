#include "kernels.h"
#ifdef __ARM_FEATURE_SVE
    #include <arm_sve.h>
#elif defined(__AVX512F__)
    #include <immintrin.h>
#elif defined(__AVX2__)
    #include <immintrin.h>
#elif defined(__AVX__)
    #include <immintrin.h>
#else
    #warning Intrinsic version not found, fall-back will be called
#endif /* __ARM_FEATURE_SVE */ 
#include "print.h"



#ifdef __ARM_FEATURE_SVE
void spmv_a64fx_sell32(densemat* b, sparsemat* mat, densemat* x)
{
    const uint64_t C_value = 32;
    if(mat->C != C_value)
    {
        printf("Wrong kernel this is SELL-%d\n", (int)C_value);
    }
    else
    {
        INFO_PRINT("SVE: SELL-%d kernel\n", C_value);
    }

    if(8 != svcntd())
    {
        printf("Wrong kernel this is 512-bit SIMD wide kernel\n");
    }

#pragma omp parallel for schedule(static)
    for(int chunk=0; chunk<mat->nchunks; ++chunk)
    {
        uint64_t idx = 0;
        uint64_t base = mat->chunkPtr[chunk];
        uint64_t n = mat->chunkLen[chunk]*C_value;

        svfloat64_t tmp0; tmp0 = svadd_z(svpfalse(),tmp0, tmp0);
        svfloat64_t tmp1; tmp1 = svadd_z(svpfalse(),tmp1, tmp1);
        svfloat64_t tmp2; tmp2 = svadd_z(svpfalse(),tmp2, tmp2);
        svfloat64_t tmp3; tmp3 = svadd_z(svpfalse(),tmp3, tmp3);
        svbool_t pg = svwhilelt_b64(idx, n);
        double *base_val0 =  &(mat->valSellC[base+0*svcntd()]);
        double *base_val1 =  &(mat->valSellC[base+1*svcntd()]);
        double *base_val2 =  &(mat->valSellC[base+2*svcntd()]);
        double *base_val3 =  &(mat->valSellC[base+3*svcntd()]);
        int *base_col0 =  &(mat->colSellC[base+0*svcntd()]);
        int *base_col1 =  &(mat->colSellC[base+1*svcntd()]);
        int *base_col2 =  &(mat->colSellC[base+2*svcntd()]);
        int *base_col3 =  &(mat->colSellC[base+3*svcntd()]);

        do
        {
            svfloat64_t mat_val0 = svld1(pg, base_val0+idx);
            svfloat64_t mat_val1 = svld1(pg, base_val1+idx);
            svfloat64_t mat_val2 = svld1(pg, base_val2+idx);
            svfloat64_t mat_val3 = svld1(pg, base_val3+idx);
            svuint64_t mat_col0 = svld1sw_u64(pg, base_col0+idx);
            svuint64_t mat_col1 = svld1sw_u64(pg, base_col1+idx);
            svuint64_t mat_col2 = svld1sw_u64(pg, base_col2+idx);
            svuint64_t mat_col3 = svld1sw_u64(pg, base_col3+idx);
            svfloat64_t x_val0 = svld1_gather_index(pg, x->val, mat_col0);
            svfloat64_t x_val1 = svld1_gather_index(pg, x->val, mat_col1);
            svfloat64_t x_val2 = svld1_gather_index(pg, x->val, mat_col2);
            svfloat64_t x_val3 = svld1_gather_index(pg, x->val, mat_col3);
            tmp0 = svmla_m(pg, tmp0, mat_val0, x_val0);
            tmp1 = svmla_m(pg, tmp1, mat_val1, x_val1);
            tmp2 = svmla_m(pg, tmp2, mat_val2, x_val2);
            tmp3 = svmla_m(pg, tmp3, mat_val3, x_val3);
            idx += 4*svcntd();
            pg = svwhilelt_b64(idx, n);
        } while(svptest_any(svptrue_b64(), pg)); 

        svst1(svptrue_b64(), &(b->val[C_value*chunk+0*svcntd()]), tmp0);
        svst1(svptrue_b64(), &(b->val[C_value*chunk+1*svcntd()]), tmp1);
        svst1(svptrue_b64(), &(b->val[C_value*chunk+2*svcntd()]), tmp2);
        svst1(svptrue_b64(), &(b->val[C_value*chunk+3*svcntd()]), tmp3);
    }
}
#elif defined(__AVX512F__)
void spmv_a64fx_sell32(densemat* b, sparsemat* mat, densemat* x)
{
    const uint64_t C_value = 32;
    if(mat->C != C_value)
    {
        ERROR_PRINT("Wrong kernel this is SELL-%d\n", (int)C_value);
    }
    else
    {
        INFO_PRINT("AVX512: SELL-%d kernel\n", C_value);
    }

#pragma omp parallel
    {
        double *rval = x->val;
        double *mval = mat->valSellC;
        int *col = mat->colSellC;
        __m512d val;
        __m512d rhs = _mm512_setzero_pd();

#pragma omp for schedule(static)
        for(int chunk=0; chunk<mat->nchunks; ++chunk)
        {
            __m512d tmp0 = _mm512_setzero_pd();
            __m512d tmp1 = _mm512_setzero_pd();
            __m512d tmp2 = _mm512_setzero_pd();
            __m512d tmp3 = _mm512_setzero_pd();
            int offs = mat->chunkPtr[chunk];
            __m256i idx;
#pragma GCC unroll 1 //no unrolling
            for(int j=0; j<mat->chunkLen[chunk]; ++j)
            {

                    val = _mm512_load_pd(&mval[offs]);idx = _mm256_load_si256((__m256i *)(&col[offs]));rhs = _mm512_i32gather_pd(idx,rval,8); tmp0 = _mm512_fmadd_pd(val,rhs,tmp0); offs+=8;
                    val = _mm512_load_pd(&mval[offs]);idx = _mm256_load_si256((__m256i *)(&col[offs]));rhs = _mm512_i32gather_pd(idx,rval,8); tmp1 = _mm512_fmadd_pd(val,rhs,tmp1); offs+=8;
                    val = _mm512_load_pd(&mval[offs]);idx = _mm256_load_si256((__m256i *)(&col[offs]));rhs = _mm512_i32gather_pd(idx,rval,8); tmp2 = _mm512_fmadd_pd(val,rhs,tmp2); offs+=8;
                    val = _mm512_load_pd(&mval[offs]);idx = _mm256_load_si256((__m256i *)(&col[offs]));rhs = _mm512_i32gather_pd(idx,rval,8); tmp3 = _mm512_fmadd_pd(val,rhs,tmp3); offs+=8;
            }
            _mm512_storeu_pd(&(b->val[chunk*C_value+(8*0)]),tmp0);
            _mm512_storeu_pd(&(b->val[chunk*C_value+(8*1)]),tmp1);
            _mm512_storeu_pd(&(b->val[chunk*C_value+(8*2)]),tmp2);
            _mm512_storeu_pd(&(b->val[chunk*C_value+(8*3)]),tmp3);
        }
    }
}

#elif defined(__AVX2__)
void spmv_a64fx_sell32(densemat* b, sparsemat* mat, densemat* x)
{
    const uint64_t C_value = 32;
    if(mat->C != C_value)
    {
        ERROR_PRINT("Wrong kernel this is SELL-%d\n", (int)C_value);
    }
    else
    {
        INFO_PRINT("AVX2: SELL-%d kernel\n", C_value);
    }


#pragma omp parallel
    {
        double *mval = mat->valSellC;
        double *rval0 = x->val;
        double *rval1 = x->val;
        double *rval2 = x->val;
        double *rval3 = x->val;
        double *rval4 = x->val;
        double *rval5 = x->val;
        double *rval6 = x->val;
        double *rval7 = x->val;
        int *col = mat->colSellC;
        double *rval = x->val;

        __m256d val;
        __m256d rhs = _mm256_setzero_pd();
#pragma omp for schedule(static)
        for(int chunk=0; chunk<mat->nchunks; ++chunk)
        {

            __m256d tmp0 = _mm256_setzero_pd();
            __m256d tmp1 = _mm256_setzero_pd();
            __m256d tmp2 = _mm256_setzero_pd();
            __m256d tmp3 = _mm256_setzero_pd();
            __m256d tmp4 = _mm256_setzero_pd();
            __m256d tmp5 = _mm256_setzero_pd();
            __m256d tmp6 = _mm256_setzero_pd();
            __m256d tmp7 = _mm256_setzero_pd();
            int offs = mat->chunkPtr[chunk];

#pragma GCC unroll 1 //no unrolling
            for(int j=0; j<mat->chunkLen[chunk]; ++j)
            {
                val    = _mm256_load_pd(&mval[offs]);rhs = _mm256_i32gather_pd(rval,_mm_load_si128((__m128i *)(col+offs)),8); tmp0   = _mm256_fmadd_pd(val,rhs,tmp0);offs+=4;
                val    = _mm256_load_pd(&mval[offs]);rhs = _mm256_i32gather_pd(rval,_mm_load_si128((__m128i *)(col+offs)),8); tmp1   = _mm256_fmadd_pd(val,rhs,tmp1);offs+=4;
                val    = _mm256_load_pd(&mval[offs]);rhs = _mm256_i32gather_pd(rval,_mm_load_si128((__m128i *)(col+offs)),8); tmp2   = _mm256_fmadd_pd(val,rhs,tmp2);offs+=4;
                val    = _mm256_load_pd(&mval[offs]);rhs = _mm256_i32gather_pd(rval,_mm_load_si128((__m128i *)(col+offs)),8); tmp3   = _mm256_fmadd_pd(val,rhs,tmp3);offs+=4;
                val    = _mm256_load_pd(&mval[offs]);rhs = _mm256_i32gather_pd(rval,_mm_load_si128((__m128i *)(col+offs)),8); tmp4   = _mm256_fmadd_pd(val,rhs,tmp4);offs+=4;
                val    = _mm256_load_pd(&mval[offs]);rhs = _mm256_i32gather_pd(rval,_mm_load_si128((__m128i *)(col+offs)),8); tmp5   = _mm256_fmadd_pd(val,rhs,tmp5);offs+=4;
                val    = _mm256_load_pd(&mval[offs]);rhs = _mm256_i32gather_pd(rval,_mm_load_si128((__m128i *)(col+offs)),8); tmp6   = _mm256_fmadd_pd(val,rhs,tmp6);offs+=4;
                val    = _mm256_load_pd(&mval[offs]);rhs = _mm256_i32gather_pd(rval,_mm_load_si128((__m128i *)(col+offs)),8); tmp7   = _mm256_fmadd_pd(val,rhs,tmp7);offs+=4;
            }
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*0)]),tmp0);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*1)]),tmp1);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*2)]),tmp2);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*3)]),tmp3);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*4)]),tmp4);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*5)]),tmp5);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*6)]),tmp6);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*7)]),tmp7);
        }
    }
}
#elif defined(__AVX__)
void spmv_a64fx_sell32(densemat* b, sparsemat* mat, densemat* x)
{
    const uint64_t C_value = 32;
    if(mat->C != C_value)
    {
        ERROR_PRINT("Wrong kernel this is SELL-%d\n", (int)C_value);
    }
    else
    {
        INFO_PRINT("AVX: SELL-%d kernel\n", C_value);
    }


#pragma omp parallel
    {
        double *mval = mat->valSellC;
        double *rval0 = x->val;
        double *rval1 = x->val;
        double *rval2 = x->val;
        double *rval3 = x->val;
        double *rval4 = x->val;
        double *rval5 = x->val;
        double *rval6 = x->val;
        double *rval7 = x->val;
        int *col = mat->colSellC;
        double *rval = x->val;

        __m256d val;
        __m256d rhs = _mm256_setzero_pd();
#pragma omp for schedule(static)
        for(int chunk=0; chunk<mat->nchunks; ++chunk)
        {

            __m256d tmp0 = _mm256_setzero_pd();
            __m256d tmp1 = _mm256_setzero_pd();
            __m256d tmp2 = _mm256_setzero_pd();
            __m256d tmp3 = _mm256_setzero_pd();
            __m256d tmp4 = _mm256_setzero_pd();
            __m256d tmp5 = _mm256_setzero_pd();
            __m256d tmp6 = _mm256_setzero_pd();
            __m256d tmp7 = _mm256_setzero_pd();
            int offs = mat->chunkPtr[chunk];

#pragma GCC unroll 1 //no unrolling
            for(int j=0; j<mat->chunkLen[chunk]; ++j)
            {
                    val = _mm256_load_pd(&mval[offs]); rhs = _mm256_setr_pd(rval[col[offs]],rval[col[offs+1]],rval[col[offs+2]],rval[col[offs+3]]); tmp0 = _mm256_add_pd(tmp0,_mm256_mul_pd(val,rhs));offs+=4;
                    val = _mm256_load_pd(&mval[offs]); rhs = _mm256_setr_pd(rval[col[offs]],rval[col[offs+1]],rval[col[offs+2]],rval[col[offs+3]]); tmp1 = _mm256_add_pd(tmp1,_mm256_mul_pd(val,rhs));offs+=4;
                    val = _mm256_load_pd(&mval[offs]); rhs = _mm256_setr_pd(rval[col[offs]],rval[col[offs+1]],rval[col[offs+2]],rval[col[offs+3]]); tmp2 = _mm256_add_pd(tmp2,_mm256_mul_pd(val,rhs));offs+=4;
                    val = _mm256_load_pd(&mval[offs]); rhs = _mm256_setr_pd(rval[col[offs]],rval[col[offs+1]],rval[col[offs+2]],rval[col[offs+3]]); tmp3 = _mm256_add_pd(tmp3,_mm256_mul_pd(val,rhs));offs+=4;
                    val = _mm256_load_pd(&mval[offs]); rhs = _mm256_setr_pd(rval[col[offs]],rval[col[offs+1]],rval[col[offs+2]],rval[col[offs+3]]); tmp4 = _mm256_add_pd(tmp4,_mm256_mul_pd(val,rhs));offs+=4;
                    val = _mm256_load_pd(&mval[offs]); rhs = _mm256_setr_pd(rval[col[offs]],rval[col[offs+1]],rval[col[offs+2]],rval[col[offs+3]]); tmp5 = _mm256_add_pd(tmp5,_mm256_mul_pd(val,rhs));offs+=4;
                    val = _mm256_load_pd(&mval[offs]); rhs = _mm256_setr_pd(rval[col[offs]],rval[col[offs+1]],rval[col[offs+2]],rval[col[offs+3]]); tmp6 = _mm256_add_pd(tmp6,_mm256_mul_pd(val,rhs));offs+=4;
                    val = _mm256_load_pd(&mval[offs]); rhs = _mm256_setr_pd(rval[col[offs]],rval[col[offs+1]],rval[col[offs+2]],rval[col[offs+3]]); tmp7 = _mm256_add_pd(tmp7,_mm256_mul_pd(val,rhs));offs+=4;

            }
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*0)]),tmp0);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*1)]),tmp1);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*2)]),tmp2);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*3)]),tmp3);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*4)]),tmp4);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*5)]),tmp5);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*6)]),tmp6);
            _mm256_storeu_pd(&(b->val[chunk*C_value+(4*7)]),tmp7);
        }
    }
}
#else
void spmv_a64fx_sell32(densemat* b, sparsemat* mat, densemat* x)
{
    //Call fallback
    spmv(b, mat, x);
}
#endif


