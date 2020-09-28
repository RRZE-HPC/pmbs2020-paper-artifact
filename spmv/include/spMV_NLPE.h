#ifndef spMV_HPP
#define spMV_HPP

#include "sparsemat.h"
#include "densemat.h"

#ifdef USE_LIKWID
extern "C"
{
#include <likwid.h>
}
#endif

#include "macros.h"
#include "print.h"

/*****CSR_MATRIX**************************************************************/
/**
 * sparse Matrix-Vector multiplication
 * y= y + A*x
 * using the CSR Format
 * y and x musst be allocated and valid
 * if accelorators are used (openACC) data have to be preent on the device
 */
void spMV_CRS_NLPE( sparsemat* A,
        densemat *x_d,
        densemat *y_d
        )
{
    double const *val     = A->val;
    int const *colInd     = A->col;
    int const *rowPtr     = A->rowPtr;
    int const numRows     = A->nrows;
    int const numNonZeros = A->nnz;

    double *x = x_d->val;
    double *y = y_d->val;
    // loop over all rows
    BEGIN_SCC(x);
#pragma omp parallel for schedule(static)
#pragma acc parallel present(val[0:numNonZeros],            \
        colInd[0:numNonZeros],      \
        rowPtr[0:numRows+1],        \
        x[0:numRows],               \
        y[0:numRows])               \
    loop
    for (int rowID=0; rowID<numRows; ++rowID)
    {
        double tmp = 0;//y[rowID];

        // loop over all elements in row
        for (int rowEntry=rowPtr[rowID]; rowEntry<rowPtr[rowID+1]; ++rowEntry)
        {
            tmp += val[rowEntry] * x[ colInd[rowEntry] ];
        }

        y[rowID] = tmp;
    }
    END_SCC;
}



/*****SELL-C-SIGMA************************************************************/

/**
 * sparse Matrix-Vector multiplication
 * y= y + A*x
 * using the Sell-C-Sigma Format
 * y and x musst be allocated and valid
 *
 * y musst be large enough to hold values for all paddded rows!
 * x must be permutaed!
 * y will be permutaed!
 */
    template< int C>
void spMV_sellC_NLPE( sparsemat* A,
        densemat *x_d,
        densemat *y_d
        )
{
    double const * val       = A->valSellC;
    int const * chunkPtr     = A->chunkPtr;
    int const * chunkLength  = A->chunkLen;
    int const * colInd       = A->colSellC;
    int const numberOfChunks = A->nchunks;
    int const chunkSize      = C;
    int const paddedRows     = A->nchunks*C;
    int const capasety       = A->nnzSellC;

    double *x = x_d->val;
    double *y = y_d->val;

    double tmp[chunkSize];


    BEGIN_SCC(x);
#pragma omp parallel for schedule(static) private(tmp)
#pragma acc parallel present(val[0 : capasety],                     \
        colInd[0 : capasety],                  \
        chunkPtr[0 : numberOfChunks],          \
        chunkLength[0 : numberOfChunks],       \
        x[0 : paddedRows],                     \
        y[0 : paddedRows])                     \
    create(tmp) private(tmp)                       \
    vector_length(32)                              \
    loop
    // loop over all chunks
    for (int chunk=0; chunk < numberOfChunks; ++chunk)
    {
        int chunkOffset = chunkPtr[chunk];

        // fill tempory vector with values from y
        for (int cRow=0        ,   rowID=chunk*chunkSize;
                cRow<chunkSize;
                ++cRow          , ++rowID
            )
        {
            tmp[cRow] = 0;//y[rowID];
        }

        // loop over all row elements in chunk
        for (int rowEntry=0; rowEntry<chunkLength[chunk]; ++rowEntry)
        {
            // (auto) vectorised loop over all rows in chunk
#pragma omp simd simdlen(8)
#pragma acc loop vector
            for (int cRow=0; cRow<chunkSize; ++cRow)
            {
                tmp[cRow] += val      [chunkOffset + rowEntry*chunkSize + cRow]
                    * x[ colInd[chunkOffset + rowEntry*chunkSize + cRow] ];
            }
        }

        // write back result of y = alpha Ax + beta y
#pragma acc loop vector
        //TODO brauch ich hier das vector und warum nihct oben?
        for (int cRow=0        , rowID=chunk*chunkSize;
                cRow<chunkSize;
                ++cRow          , ++rowID
            )
        {
            y[rowID] = tmp[cRow];
        }
    }
    END_SCC;
}

/*wrapper function for dynamic dispatshing*/
void spMV_NLPE( densemat *y,
        sparsemat *A,
        densemat *x
        )
{
    int C = A->C;

    if (1 == C)
        return spMV_CRS_NLPE(A,x,y);
    else if (4 == C)
        return spMV_sellC_NLPE<4>(A,x,y);
    else if (8 == C)
        return spMV_sellC_NLPE<8>(A,x,y);
    else if (16 == C)
        return spMV_sellC_NLPE<16>(A,x,y);
    else if (24 == C)
        return spMV_sellC_NLPE<24>(A,x,y);
    else if (32 == C)
        return spMV_sellC_NLPE<32>(A,x,y);
    else if (64 == C)
        return spMV_sellC_NLPE<64>(A,x,y);
    else if (128 == C)
        return spMV_sellC_NLPE<128>(A,x,y);

#ifdef SET_C
    else if (SET_C == C)
        return spMV_sellC_NLPE<SET_C>(A,x,y);
#endif
    else
    {
        WARNING_PRINT("No specialized kernel found for NLPE impl\n");
        return spMV_sellC_NLPE<1>(A,x,y);
    }
}


#endif
