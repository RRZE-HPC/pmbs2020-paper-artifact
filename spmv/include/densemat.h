#ifndef RACE_DENSEMAT_H
#define RACE_DENSEMAT_H

#include <functional>

struct densemat
{
    int nrows;
    int ncols;
    double *val;

    void setVal(double value);
    void setRand();
//    void setFn(std::function<double(int)> fn);
    void setFn(std::function<double(void)> fn);

    densemat(int nrows, int ncols);
    ~densemat();

};

bool checkEqual(const densemat* lhs, const densemat* rhs, double tol=1e-4);

#endif
