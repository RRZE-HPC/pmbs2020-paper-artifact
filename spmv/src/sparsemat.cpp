#include "sparsemat.h"
#include "mmio.h"
#include "stdlib.h"
#include <omp.h>
#include <vector>
#include <algorithm>
#ifdef USE_SPMP
#include "CSR.hpp"
#include "reordering/BFSBipartite.hpp"
#endif
#include "kernels.h"
#include "timer.h"

sparsemat::sparsemat():nrows(0), nnz(0), val(NULL), rowPtr(NULL), col(NULL), nnz_symm(0), rowPtr_symm(NULL), col_symm(NULL), val_symm(NULL), chunkLen(NULL), chunkPtr(NULL), colSellC(NULL), valSellC(NULL), unrollFac(1), C(1)
{

#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    rcmPerm = NULL;
    rcmInvPerm = NULL;
    nnzPerRow = NULL;
}

sparsemat::~sparsemat()
{
    if(val)
        delete[] val;

    if(rowPtr)
        delete[] rowPtr;

    if(col)
        delete[] col;

    if(val_symm)
        delete[] val_symm;

    if(rowPtr_symm)
        delete[] rowPtr_symm;

    if(col_symm)
        delete[] col_symm;

    if(chunkLen)
        delete[] chunkLen;

    if(chunkPtr)
        delete[] chunkPtr;

    if(colSellC)
        delete[] colSellC;

    if(valSellC)
        delete[] valSellC;

    if(nnzPerRow)
    {
        delete[] nnzPerRow;
    }
}

void sparsemat::generateHPCG(int nx, int ny, int nz)
{

    int numberOfNonzerosPerRow = 27; // We are approximating a 27-point finite element/volume/difference 3D stencil
    int numberOfRows = nx*ny*nz;

    // Allocate arrays that are of length localNumberOfRows
    int * nonzerosInRow = new int[numberOfRows];
    int ** mtxInd = new int*[numberOfRows];
    double ** matrixValues = new double*[numberOfRows];

    int *col_ = new int[numberOfNonzerosPerRow*numberOfRows];
    double *val_ = new double[sizeof(double)*numberOfNonzerosPerRow*numberOfRows];

    int* boundaryRows = new int[nx*ny*nz - (nx-2)*(ny-2)*(nz-2)];

    if ( col_ == NULL || val_ == NULL || nonzerosInRow == NULL || boundaryRows == NULL
            || mtxInd == NULL || matrixValues == NULL)
    {
        return;
    }

    int numOfBoundaryRows = 0;
#pragma omp parallel reduction(+:numOfBoundaryRows)
    {
#pragma omp for nowait
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                boundaryRows[y*nx + x] = y*nx + x;
                numOfBoundaryRows++;
            }
        }

#pragma omp for nowait
        for (int z = 1; z < nz - 1; z++) {
            for (int x = 0; x < nx; x++) {
                boundaryRows[ny*nx + 2*(z-1)*(nx+ny-2) + x ] = z*ny*nx + x;
                numOfBoundaryRows++;
            }
            for (int y = 1; y < ny - 1; y++) {
                boundaryRows[ny*nx + 2*(z-1)*(nx+ny-2) + nx + 2*(y-1)] = (z*ny + y)*nx;
                numOfBoundaryRows++;
                boundaryRows[ny*nx + 2*(z-1)*(nx+ny-2) + nx + 2*(y-1)+1] = (z*ny + y)*nx + nx - 1;
                numOfBoundaryRows++;
            }
            for (int x = 0; x < nx; x++) {
                boundaryRows[ny*nx + 2*(z-1)*(nx+ny-2) + nx + 2*(ny-2) + x] = (z*ny + (ny - 1))*nx + x;
                numOfBoundaryRows++;
            }
        }

#pragma omp for nowait
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                boundaryRows[ny*nx + 2*(nz-2)*(nx+ny-2) + y*nx + x] = ((nz - 1)*ny + y)*nx + x;
                numOfBoundaryRows++;
            }
        }
    }

    int numberOfNonzeros = 0;



#pragma omp parallel reduction(+:numberOfNonzeros)
    {
        int ithr = omp_get_thread_num();
        int nthr = omp_get_num_threads();

        int works = (nz - 2)*(ny - 2);
        int begin = ((ithr  )*works)/nthr;
        int end   = ((ithr+1)*works)/nthr;
        for (int i = begin; i < end; i++)
        {
            int iz = i/(ny - 2) + 1;
            int iy = i%(ny - 2) + 1;

            for (int ix=1; ix<nx-1; ix++)
            {
                int currentLocalRow = iz*nx*ny+iy*nx+ix;
                mtxInd[currentLocalRow]      = col_ + currentLocalRow*numberOfNonzerosPerRow;
                matrixValues[currentLocalRow] = val_ + currentLocalRow*numberOfNonzerosPerRow;
                char numberOfNonzerosInRow = 0;
                double * currentValuePointer = matrixValues[currentLocalRow]; // Pointer to current value in current row
                int  * currentIndexPointerL = mtxInd[currentLocalRow]; // Pointer to current index in current row
                for (int sz=-1; sz<=1; sz++) {

                    *(currentValuePointer + 0) = -1.0;
                    *(currentValuePointer + 1) = -1.0;
                    *(currentValuePointer + 2) = -1.0;
                    *(currentValuePointer + 3) = -1.0;
                    *(currentValuePointer + 4) = -1.0;
                    *(currentValuePointer + 5) = -1.0;
                    *(currentValuePointer + 6) = -1.0;
                    *(currentValuePointer + 7) = -1.0;
                    *(currentValuePointer + 8) = -1.0;

                    int offset = currentLocalRow + sz*ny*nx;
                    *(currentIndexPointerL + 0) = offset - nx - 1;
                    *(currentIndexPointerL + 1) = offset - nx;
                    *(currentIndexPointerL + 2) = offset - nx + 1;
                    *(currentIndexPointerL + 3) = offset - 1;
                    *(currentIndexPointerL + 4) = offset;
                    *(currentIndexPointerL + 5) = offset + 1;
                    *(currentIndexPointerL + 6) = offset + nx - 1;
                    *(currentIndexPointerL + 7) = offset + nx;
                    *(currentIndexPointerL + 8) = offset + nx + 1;

                    currentValuePointer  += 9;
                    currentIndexPointerL += 9;
                } // end sz loop
                *(currentValuePointer - 14) = 26.0;
                numberOfNonzerosInRow += 27;
                nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow;
                numberOfNonzeros += numberOfNonzerosInRow; // Protect this with an atomic
            } // end ix loop
        }

#pragma omp for
        for (int i = 0; i < numOfBoundaryRows; i++) {
            int currentLocalRow = boundaryRows[i];

            int iz = currentLocalRow/(ny*nx);
            int iy = currentLocalRow/nx%ny;
            int ix = currentLocalRow%nx;

            int sz_begin = std::max<int>(-1, -iz);
            int sz_end = std::min<int>(1, nz - iz - 1);

            int sy_begin = std::max<int>(-1, -iy);
            int sy_end = std::min<int>(1, ny - iy - 1);

            int sx_begin = std::max<int>(-1, -ix);
            int sx_end = std::min<int>(1, nx - ix - 1);


            mtxInd[currentLocalRow]      = col_ + currentLocalRow*numberOfNonzerosPerRow;
            matrixValues[currentLocalRow] = val_ + currentLocalRow*numberOfNonzerosPerRow;
            char numberOfNonzerosInRow = 0;
            double * currentValuePointer = matrixValues[currentLocalRow]; // Pointer to current value in current row
            int  * currentIndexPointerL = mtxInd[currentLocalRow];
            for (int sz=sz_begin; sz<=sz_end; sz++) {
                for (int sy=sy_begin; sy<=sy_end; sy++) {
                    for (int sx=sx_begin; sx<=sx_end; sx++) {
                        int    col = currentLocalRow + sz*nx*ny+sy*nx+sx;
                        if (col==currentLocalRow) {
                            *currentValuePointer++ = 26.0;
                        } else {
                            *currentValuePointer++ = -1.0;
                        }
                        *currentIndexPointerL++ = col;
                        numberOfNonzerosInRow++;
                    } // end sx loop
                } // end sy loop
            } // end sz loop
            nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow;
            numberOfNonzeros += numberOfNonzerosInRow; // Protect this with an atomic
        }
    }

    rowPtr = new int[numberOfRows+1];
    col = new int[numberOfNonzeros];
    val = new double[numberOfNonzeros];

    rowPtr[0] = 0;
#pragma omp parallel for
    for(int i=0; i<numberOfRows; ++i)
    {
        rowPtr[i+1] = nonzerosInRow[i];
    }

    for(int i=0; i<numberOfRows; ++i)
    {
        rowPtr[i+1] += rowPtr[i];
    }

#pragma omp parallel for
    for(int i=0; i<numberOfRows; ++i)
    {
        int k = rowPtr[i];
        for(int j=0; j<nonzerosInRow[i];++j)
        {
            col[k] = (mtxInd[i])[j];
            val[k] = (matrixValues[i])[j];
            ++k;
        }
    }

    nrows = numberOfRows;
    nnz = numberOfNonzeros;

    delete [] col_;
    delete [] val_;
    delete [] mtxInd;
    delete [] matrixValues;
    nnzPerRow = nonzerosInRow;
//    free(nonzerosInRow);

    return;
}


bool sparsemat::readFile(char* filename)
{

    MM_typecode matcode;

    FILE *f;

    if ((f = fopen(filename, "r")) == NULL)
        return -1;


    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner ");
        printf(" in file [%s]\n", filename);
        return -1;
    }

    fclose(f);

    bool compatible_flag = (mm_is_sparse(matcode) && (mm_is_real(matcode)||mm_is_pattern(matcode))) && (mm_is_symmetric(matcode) || mm_is_general(matcode));
    bool symm_flag = mm_is_symmetric(matcode);
    bool pattern_flag = mm_is_pattern(matcode);

    if(!compatible_flag)
    {
        printf("The matrix market file provided is not supported.\n Reason :\n");
        if(!mm_is_sparse(matcode))
        {
            printf(" * matrix has to be sparse\n");
        }

        if(!mm_is_real(matcode) && !(mm_is_pattern(matcode)))
        {
            printf(" * matrix has to be real or pattern\n");
        }

        if(!mm_is_symmetric(matcode) && !mm_is_general(matcode))
        {
            printf(" * matrix has to be either general or symmetric\n");
        }

        exit(0);
    }

    int ncols;
    int *row;
    int *col_unsorted;
    double *val_unsorted;

    if(mm_read_unsymmetric_sparse(filename, &nrows, &ncols, &nnz, &val_unsorted, &row, &col_unsorted) < 0)
    {
        printf("Error in file reading\n");
        return false;
    }
    if(nrows != ncols)
    {
        printf("Currently only Symmetric matrices are supported\n");
        return false;
    }

    //If matrix market file is symmetric; create a general one out of it
    if(symm_flag)
    {
        printf("Creating a general matrix out of a symmetric one\n");

        int ctr = 0;

        //this is needed since diagonals might be missing in some cases
        for(int idx=0; idx<nnz; ++idx)
        {
            ++ctr;
            if(row[idx]!=col_unsorted[idx])
            {
                ++ctr;
            }
        }

        int new_nnz = ctr;

        int *row_general = new int[new_nnz];
        int *col_general = new int[new_nnz];
        double *val_general = new double[new_nnz];

        int idx_gen=0;

        for(int idx=0; idx<nnz; ++idx)
        {
            row_general[idx_gen] = row[idx];
            col_general[idx_gen] = col_unsorted[idx];
            val_general[idx_gen] = val_unsorted[idx];
            ++idx_gen;

            if(row[idx] != col_unsorted[idx])
            {
                row_general[idx_gen] = col_unsorted[idx];
                col_general[idx_gen] = row[idx];
                val_general[idx_gen] = val_unsorted[idx];
                ++idx_gen;
            }
        }

        free(row);
        free(col_unsorted);
        free(val_unsorted);

        nnz = new_nnz;

        //assign right pointers for further proccesing
        row = row_general;
        col_unsorted = col_general;
        val_unsorted = val_general;
    }

    //permute the col and val according to row
    int* perm = new int[nnz];
    for(int idx=0; idx<nnz; ++idx)
    {
        perm[idx] = idx;
    }

    sort_perm(row, perm, nnz);

    col = new int[nnz];
    val = new double[nnz];

    for(int idx=0; idx<nnz; ++idx)
    {
        col[idx] = col_unsorted[perm[idx]];
        val[idx] = val_unsorted[perm[idx]];
    }

    delete[] col_unsorted;
    delete[] val_unsorted;


    rowPtr = new int[nrows+1];

    nnzPerRow = new int[nrows];
    for(int i=0; i<nrows; ++i)
    {
        nnzPerRow[i] = 0;
    }

    //count nnz per row
    for(int i=0; i<nnz; ++i)
    {
        ++nnzPerRow[row[i]];
    }

    rowPtr[0] = 0;
    for(int i=0; i<nrows; ++i)
    {
        rowPtr[i+1] = rowPtr[i]+nnzPerRow[i];
    }

    if(rowPtr[nrows] != nnz)
    {
        printf("Error in reading matrix\n");
        return false;
    }

    delete[] row;
    delete[] perm;

    printf("NUMA initing\n");
    NUMAinit();
    printf("NUMA inited\n");
    //writeFile("beforePerm.mtx");
    return true;
}


#if 0
bool sparsemat::readFile(char* filename)
{

    MM_typecode matcode;

    FILE *f;

    if ((f = fopen(filename, "r")) == NULL)
        return -1;


    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner ");
        printf(" in file [%s]\n", filename);
        return -1;
    }

    fclose(f);

    bool compatible_flag = (mm_is_sparse(matcode) && (mm_is_real(matcode)||mm_is_pattern(matcode))) && (mm_is_symmetric(matcode) || mm_is_general(matcode));
    bool symm_flag = mm_is_symmetric(matcode);
    bool pattern_flag = mm_is_pattern(matcode);

    if(!compatible_flag)
    {
        printf("The matrix market file provided is not supported.\n Reason :\n");
        if(!mm_is_sparse(matcode))
        {
            printf(" * matrix has to be sparse\n");
        }

        if(!mm_is_real(matcode) && !(mm_is_pattern(matcode)))
        {
            printf(" * matrix has to be real or pattern\n");
        }

        if(!mm_is_symmetric(matcode) && !mm_is_general(matcode))
        {
            printf(" * matrix has to be either general or symmetric\n");
        }

        exit(0);
    }

    int ncols;
    int *row;
    int *col_unsorted;
    double *val_unsorted;

    if(mm_read_unsymmetric_sparse(filename, &nrows, &ncols, &nnz, &val_unsorted, &row, &col_unsorted) < 0)
    {
        printf("Error in file reading\n");
        return false;
    }
    if(nrows != ncols)
    {
        printf("Currently only Symmetric matrices are supported\n");
        return false;
    }

    //If matrix market file is symmetric; create a general one out of it
    if(symm_flag)
    {
        printf("Creating a general matrix out of a symmetric one\n");

        int ctr = 0;

        //this is needed since diagonals might be missing in some cases
        for(int idx=0; idx<nnz; ++idx)
        {
            ++ctr;
            if(row[idx]!=col_unsorted[idx])
            {
                ++ctr;
            }
        }

        int new_nnz = ctr;

        int *row_general = new int[new_nnz];
        int *col_general = new int[new_nnz];
        double *val_general = new double[new_nnz];

        int idx_gen=0;

        for(int idx=0; idx<nnz; ++idx)
        {
            row_general[idx_gen] = row[idx];
            col_general[idx_gen] = col_unsorted[idx];
            val_general[idx_gen] = val_unsorted[idx];
            ++idx_gen;

            if(row[idx] != col_unsorted[idx])
            {
                row_general[idx_gen] = col_unsorted[idx];
                col_general[idx_gen] = row[idx];
                val_general[idx_gen] = val_unsorted[idx];
                ++idx_gen;
            }
        }

        free(row);
        free(col_unsorted);
        free(val_unsorted);

        nnz = new_nnz;

        //assign right pointers for further proccesing
        row = row_general;
        col_unsorted = col_general;
        val_unsorted = val_general;
    }

    //permute the col and val according to row
    int* perm = new int[nnz];
    for(int idx=0; idx<nnz; ++idx)
    {
        perm[idx] = idx;
    }

    sort_perm(row, perm, nnz);

    col = new int[nnz];
    val = new double[nnz];

    for(int idx=0; idx<nnz; ++idx)
    {
        col[idx] = col_unsorted[perm[idx]];
        val[idx] = val_unsorted[perm[idx]];
    }

    delete[] col_unsorted;
    delete[] val_unsorted;


    rowPtr = new int[nrows+1];

    int *nnzPerRow = new int[nrows];
    for(int i=0; i<nrows; ++i)
    {
        nnzPerRow[i] = 0;
    }

    //count nnz per row
    for(int i=0; i<nnz; ++i)
    {
        ++nnzPerRow[row[i]];
    }

    rowPtr[0] = 0;
    for(int i=0; i<nrows; ++i)
    {
        rowPtr[i+1] = rowPtr[i]+nnzPerRow[i];
    }

    if(rowPtr[nrows] != nnz)
    {
        printf("Error in reading matrix\n");
        return false;
    }

    delete[] row;

    return true;
}
#endif

//If diag is not present fill with value given
void sparsemat::ensureDiag(double diag_val)
{
    //check whether a new allocation is necessary
    int extra_nnz=0;
    std::vector<double>* val_with_diag = new std::vector<double>();
    std::vector<int>* col_with_diag = new std::vector<int>();
    std::vector<int>* rowPtr_with_diag = new std::vector<int>(rowPtr, rowPtr+nrows+1);

    for(int row=0; row<nrows; ++row)
    {
        std::vector<double> val_row, col_row;
        bool diagHit = false;
        for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
        {
            val_row.push_back(val[idx]);
            col_row.push_back(col[idx]);

            if(col[idx] == row)
            {
                diagHit = true;
            }
        }
        if(!diagHit)
        {
            val_row.push_back(diag_val);
            col_row.push_back(row);

            //sort val and col
            std::stable_sort(col_row.begin(), col_row.end());
            std::stable_sort(val_row.begin(), val_row.end(), [&](const int& a, const int& b) {return (col_row[a] < col_row[b]); });
            ++extra_nnz;
        }
        val_with_diag->insert(val_with_diag->end(), val_row.begin(), val_row.end());
        col_with_diag->insert(col_with_diag->end(), col_row.begin(), col_row.end());
        rowPtr_with_diag->at(row+1) = rowPtr_with_diag->at(row+1) + extra_nnz;
    }
    //allocate new matrix if necessary
    if(extra_nnz)
    {
        delete[] val;
        delete[] col;
        delete[] rowPtr;

        nnz += extra_nnz;
        val = new double[nnz];
        col = new int[nnz];
        rowPtr = new int[nrows+1];

        rowPtr[0] = rowPtr_with_diag->at(0);
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row)
        {
            rowPtr[row+1] = rowPtr_with_diag->at(row+1);
            for(int idx=rowPtr_with_diag->at(row); idx<rowPtr_with_diag->at(row+1); ++idx)
            {
                val[idx] = val_with_diag->at(idx);
                col[idx] = col_with_diag->at(idx);
            }
        }
        printf("Added diagonal entries\n");
    }

    delete val_with_diag;
    delete col_with_diag;
    delete rowPtr_with_diag;
}

//necessary for GS like kernels
void sparsemat::makeDiagFirst()
{
    //check whether a new allocation is necessary
    int extra_nnz=0;
    std::vector<double>* val_with_diag = new std::vector<double>();
    std::vector<int>* col_with_diag = new std::vector<int>();
    std::vector<int>* rowPtr_with_diag = new std::vector<int>(rowPtr, rowPtr+nrows+1);

    for(int row=0; row<nrows; ++row)
    {
        bool diagHit = false;
        for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
        {
            val_with_diag->push_back(val[idx]);
            col_with_diag->push_back(col[idx]);

            if(col[idx] == row)
            {
                diagHit = true;
            }
        }
        if(!diagHit)
        {
            val_with_diag->push_back(0.0);
            col_with_diag->push_back(row);
            ++extra_nnz;
        }
        rowPtr_with_diag->at(row+1) = rowPtr_with_diag->at(row+1) + extra_nnz;
    }

    //allocate new matrix if necessary
    if(extra_nnz)
    {
        delete[] val;
        delete[] col;
        delete[] rowPtr;

        nnz += extra_nnz;
        val = new double[nnz];
        col = new int[nnz];
        rowPtr = new int[nrows+1];

        rowPtr[0] = rowPtr_with_diag->at(0);
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row)
        {
            rowPtr[row+1] = rowPtr_with_diag->at(row+1);
            for(int idx=rowPtr_with_diag->at(row); idx<rowPtr_with_diag->at(row+1); ++idx)
            {
                val[idx] = val_with_diag->at(idx);
                col[idx] = col_with_diag->at(idx);
            }
        }
        printf("Explicit 0 in diagonal entries added\n");
    }

    delete val_with_diag;
    delete col_with_diag;
    delete rowPtr_with_diag;

#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        bool diag_hit = false;

        double* newVal = new double[rowPtr[row+1]-rowPtr[row]];
        int* newCol = new int[rowPtr[row+1]-rowPtr[row]];
        for(int idx=rowPtr[row], locIdx=0; idx<rowPtr[row+1]; ++idx, ++locIdx)
        {
            //shift all elements+1 until diag entry
            if(col[idx] == row)
            {
                newVal[0] = val[idx];
                newCol[0] = col[idx];
                diag_hit = true;
            }
            else if(!diag_hit)
            {
                newVal[locIdx+1] = val[idx];
                newCol[locIdx+1] = col[idx];
            }
            else
            {
                newVal[locIdx] = val[idx];
                newCol[locIdx] = col[idx];
            }
        }
        //assign new Val
        for(int idx = rowPtr[row], locIdx=0; idx<rowPtr[row+1]; ++idx, ++locIdx)
        {
            val[idx] = newVal[locIdx];
            col[idx] = newCol[locIdx];
        }

        delete[] newVal;
        delete[] newCol;
    }

}

//write matrix market file
bool sparsemat::writeFile(char* filename)
{
    int* row_1_based = new int[nnz];
    int* col_1_based = new int[nnz];

    //create row indices
    for(int i=0; i<nrows; ++i)
    {
        for(int idx=rowPtr[i]; idx<rowPtr[i+1]; ++idx)
        {
            row_1_based[idx]=i+1;
            col_1_based[idx]=col[idx]+1;
        }
    }

    mm_write_mtx_crd(filename, nrows, nrows, nnz, row_1_based, col_1_based, val, "MCRG");

    delete[] row_1_based;
    delete[] col_1_based;
}

bool sparsemat::computeSymmData()
{
    //compute only if previously not computed
    if(nnz_symm == 0)
    {
        /* Here we compute symmetric data of matrix
         * which is used if necessary; upper symmetric
         * portion is stored*/

        nnz_symm = 0;
        rowPtr_symm = new int[nrows+1];
        rowPtr_symm[0] = 0;

        //NUMA init
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row)
        {
            rowPtr_symm[row+1] = 0;
        }

        //count non-zeros in upper-symm
        for(int row=0; row<nrows; ++row) {
            for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx) {
                if(col[idx]>=row) {
                    ++nnz_symm;
                }
                rowPtr_symm[row+1] = nnz_symm;
            }
        }

        col_symm = new int[nnz_symm];
        val_symm = new double[nnz_symm];

        //With NUMA init
#pragma omp parallel for schedule(static)
        for(int row=0; row<nrows; ++row) {
            int idx_symm = rowPtr_symm[row];
            for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx) {
                if(col[idx]>=row) {
                    val_symm[idx_symm] = val[idx];
                    col_symm[idx_symm] = col[idx];
                    ++idx_symm;
                }
            }
        }
    }
    return true;
}

void sparsemat::doRCM()
{
#ifdef USE_SPMP
    int orig_threads = 1;
    printf("Doing RCM permutation\n");
#pragma omp parallel
    {
        orig_threads = omp_get_num_threads();
    }
    omp_set_num_threads(1);

    SpMP::CSR *csr = NULL;
    csr = new SpMP::CSR(nrows, nrows, rowPtr, col, val);
    //   rcmPerm = new int[nrows];
    //    rcmInvPerm = new int[nrows];
    if(csr->isSymmetric(false,false))
    {
        rcmPerm = new int[nrows];
        rcmInvPerm = new int[nrows];
        csr->getRCMPermutation(rcmInvPerm, rcmPerm);
    }
    else
    {
        printf("Matrix not symmetric RCM cannot be done\n");
    }
    omp_set_num_threads(orig_threads);
    delete csr;
#else
    printf("RCM not enabled, please install SpMP and set ENABLE_SPMP in config.mk\n");
#endif
}

void sparsemat::doRCMPermute()
{
    doRCM();
    if(rcmPerm)
    {
        permute(rcmPerm, rcmInvPerm);
    }
    rcmPerm = NULL;
    rcmInvPerm = NULL;
}

//symmetrically permute
void sparsemat::permute(int *perm, int*  invPerm)
{
    double* newVal = new double[nnz];
    int* newRowPtr = new int[nrows+1];
    int* newCol = new int[nnz];

    newRowPtr[0] = 0;

    //NUMA init
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        newRowPtr[row+1] = 0;
    }

    //first find newRowPtr; therefore we can do proper NUMA init
    int permIdx=0;
    printf("nrows = %d\n", nrows);
    for(int row=0; row<nrows; ++row)
    {
        //row permutation
        int permRow = perm[row];
        nnzPerRow[row] = (rowPtr[permRow+1]-rowPtr[permRow]);
        for(int idx=rowPtr[permRow]; idx<rowPtr[permRow+1]; ++idx)
        {
            ++permIdx;
        }
        newRowPtr[row+1] = permIdx;
    }


    //with NUMA init
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        //row permutation
        int permRow = perm[row];
        for(int permIdx=newRowPtr[row],idx=rowPtr[permRow]; permIdx<newRowPtr[row+1]; ++idx,++permIdx)
        {
            //permute column-wise also
            newVal[permIdx] = val[idx];
            newCol[permIdx] = invPerm[col[idx]];
        }
    }

    //free old permutations
    delete[] val;
    delete[] rowPtr;
    delete[] col;

    val = newVal;
    rowPtr = newRowPtr;
    col = newCol;
}

void sparsemat::constructSellCSigma(int chunkHeight, int sigma, int pad)
{
    C = chunkHeight;
    P = pad;

    int nSigmaChunks = (int)(nrows/(double)sigma);
    if(sigma > 1)
    {
        int *sigmaPerm = new int[nrows];
        for(int i=0; i<nrows; ++i)
        {
            sigmaPerm[i] = i;
        }

        for(int sigmaChunk=0; sigmaChunk<nSigmaChunks; ++sigmaChunk)
        {
            int *perm_begin = &(sigmaPerm[sigmaChunk*sigma]);
            sort_perm(nnzPerRow, perm_begin, sigma);
        }

        int restSigmaChunk = nrows%sigma;
        if(restSigmaChunk > C)
        {
            int *perm_begin = &(sigmaPerm[nSigmaChunks*sigma]);
            sort_perm(nnzPerRow, perm_begin, restSigmaChunk);
        }

        int *sigmaInvPerm = new int[nrows];

        for(int i=0; i<nrows; ++i)
        {
            sigmaInvPerm[sigmaPerm[i]] = i;
        }

        permute(sigmaPerm, sigmaInvPerm);

        delete[] sigmaPerm;
        delete[] sigmaInvPerm;
    }

    nchunks = (int)(nrows/(double)C);
    if(nrows%C > 0)
    {
        nchunks += 1;
    }

    chunkLen = new int[nchunks];
    chunkPtr = new int[nchunks+1];

#pragma omp parallel for schedule(static)
    for(int i=0; i<nchunks; ++i)
    {
        chunkLen[i] = 0;
        chunkPtr[i] = 0;
    }

    nnzSellC = 0;
    //find chunkLen
    for(int chunk=0; chunk<nchunks; ++chunk)
    {
        int maxRowLen = 0;
        for(int rowInChunk=0; rowInChunk<C; ++rowInChunk)
        {
            int row = chunk*C + rowInChunk;
            if(row<nrows)
            {
                maxRowLen = std::max(maxRowLen, rowPtr[row+1]-rowPtr[row]);
            }
        }
        //pad it to be multiple of P
        if((maxRowLen%P) != 0)
        {
            maxRowLen = ((int)(maxRowLen/(double)P)+1)*P;
        }
        chunkLen[chunk] = maxRowLen;
        nnzSellC += maxRowLen*C;
    }

    colSellC = new int[nnzSellC];
    valSellC = new double[nnzSellC];


#pragma omp parallel for schedule(static)
    for(int i=0; i<=(nchunks); ++i)
    {
        chunkPtr[i] = 0;
    }

    for(int i=0; i<(nchunks); ++i)
    {
        chunkPtr[i+1] = chunkPtr[i] + C*chunkLen[i];
    }


#pragma omp parallel for schedule(static)
    for(int chunk=0; chunk<nchunks; ++chunk)
    {
        for(int rowInChunk=0; rowInChunk<C; ++rowInChunk)
        {
            for(int idx=0; idx<chunkLen[chunk]; ++idx)
            {
                if(C == nrows)
                {
                    colSellC[chunkPtr[chunk]+idx*C+rowInChunk] = chunk*C + rowInChunk;//ELLPACK of Fujitsu needs it this way (the rowIndex)
                }
                else
                {
                    colSellC[chunkPtr[chunk]+idx*C+rowInChunk] = 0;
                }
                valSellC[chunkPtr[chunk]+idx*C+rowInChunk] = 0;
            }
        }
    }


    for(int chunk=0; chunk<nchunks; ++chunk)
    {
        for(int rowInChunk=0; rowInChunk<C; ++rowInChunk)
        {
            int row = chunk*C + rowInChunk;
            if(row<nrows)
            {
                for(int idx=rowPtr[row],j=0; idx<rowPtr[row+1]; ++idx,++j)
                {
                    valSellC[chunkPtr[chunk]+j*C+rowInChunk] = val[idx];
                    colSellC[chunkPtr[chunk]+j*C+rowInChunk] = col[idx];
                }
            }
        }
    }

    std::vector<double> strideAvg(nchunks*C, 0);
    double strideAvg_total = 0;
    for(int chunk=0; chunk<nchunks; ++chunk)
    {
        for(int rowInChunk=0; rowInChunk<C; ++rowInChunk)
        {
            for(int idx=1; idx<chunkLen[chunk]; ++idx)
            {
                strideAvg[chunk*C+rowInChunk] += std::abs(colSellC[chunkPtr[chunk]+idx*C+rowInChunk] - colSellC[chunkPtr[chunk]+(idx-1)*C+rowInChunk]);
            }

            strideAvg[chunk*C+rowInChunk] = strideAvg[chunk*C+rowInChunk]/(double)chunkLen[chunk];
            strideAvg_total += strideAvg[chunk*C+rowInChunk];
        }
    }
    strideAvg_total = strideAvg_total/((double)nchunks*C);

    printf("Average stride length = %f\n", strideAvg_total);
}


//symmetrically permute
void sparsemat::NUMAinit()
{
    double* newVal = new double[nnz];
    int* newCol = new int[nnz];
    int* newRowPtr = new int[nrows+1];

    /*
       double *newVal = (double*) allocate(1024, sizeof(double)*nnz);
       int *newCol = (int*) allocate(1024, sizeof(int)*nnz);
       int *newRowPtr = (int*) allocate(1024, sizeof(int)*(nrows+1));
       */

    //NUMA init
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows+1; ++row)
    {
        newRowPtr[row] = rowPtr[row];
    }
#pragma omp parallel for schedule(static)
    for(int row=0; row<nrows; ++row)
    {
        for(int idx=newRowPtr[row]; idx<newRowPtr[row+1]; ++idx)
        {
            //newVal[idx] = val[idx];
            newCol[idx] = col[idx];
            newVal[idx] = val[idx];
        }
    }


    //free old _perm_utations
    delete[] val;
    delete[] rowPtr;
    delete[] col;

    val = newVal;
    rowPtr = newRowPtr;
    col = newCol;
}

