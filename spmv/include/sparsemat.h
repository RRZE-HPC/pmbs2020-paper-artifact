#ifndef RACE_SPARSEMAT_H
#define RACE_SPARSEMAT_H

#include <algorithm>
#include <iterator>
#include <vector>

template <typename T> void sort_perm(T *arr, int *perm, int len, bool rev=false)
{
    if(rev == false) {
        std::stable_sort(perm+0, perm+len, [&](const int& a, const int& b) {return (arr[a] < arr[b]);});
    } else {
        std::stable_sort(perm+0, perm+len, [&](const int& a, const int& b) {return (arr[a] > arr[b]); });
    }
}

//for debugging
template <typename T> void sort_perm_v(T *arr, int *perm, int len, bool rev=false)
{
    if(rev == false) {
        std::stable_sort(perm+0, perm+len, [&](const int& a, const int& b) {printf("comparing arr[%d] = %d, arr[%d] = %d\n", a, arr[a], b, arr[b]); return (arr[a] < arr[b]);});
    } else {
        std::stable_sort(perm+0, perm+len, [&](const int& a, const int& b) {printf("comparing arr[%d] = %d, arr[%d] = %d\n", a, arr[a], b, arr[b]); return (arr[a] > arr[b]); });
    }
}


struct sparsemat
{
    int nrows, nnz;
    //interface to coloring engine
    int *rowPtr, *col, *nnzPerRow;
    double *val;
    int *rcmPerm, *rcmInvPerm;
    int C;
    int P;
    int nchunks, nnzSellC;
    int *chunkLen;
    int *chunkPtr;
    int *colSellC;
    double *valSellC;
    int unrollFac; //for kernel, just a work-around
    int nthreads;

    bool readFile(char* filename);
    bool writeFile(char* filename);
    void ensureDiag(double diag_val=0);
    void makeDiagFirst();
    void generateHPCG(int nx, int ny, int nz);
    void doRCM();
    void doRCMPermute();
    void permute(int* perm, int* invPerm);
    void NUMAinit();

    /* For symmetric computations */
    int nnz_symm;
    int *rowPtr_symm, *col_symm;
    double *val_symm;
    bool computeSymmData();

    void constructSellCSigma(int C, int sigma, int P=1);
    sparsemat();
    ~sparsemat();

};

#endif
