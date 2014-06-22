#ifndef NULLSOLVE_INL
#define NULLSOLVE_INL

#include "NuLLSolve.h"

namespace NuLLMath
{
    template <typename T> T dotProduct(const Vector<T>& x, const Vector<T>& y);
    template <typename T> T euclideanNorm(const Vector<T>& x);
}

//for solving linear systems of equations
//use forward substitution with decomposed matrix first, then backwardsubstitution
//applicable to cholesky matrix and LU-decomposed matrix
template<typename T> void NuLLSolve::forwardElimination(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res)
{
    T sum;
    size_t dim = A.height();
    res = b;

    //first entry
    res[0] = res[0] / A(0,0);

    //remaining entries
    for(uint i=1; i<dim; ++i)
    {
        //res[i-1] = res[i-1];
        sum = 0;
        for(uint k=0; k<=i-1; ++k)
            sum += A(i,k) * res[k];

        res[i] -= sum;
        res[i] /= A(i,i);
    }
}

//for solving linear systems of equations
//use forward substitution with decomposed matrix first, then backwardsubstitution
//applicable on cholesky matrix and LU-decomposed matrix
template<typename T> void NuLLSolve::backwardSubstitution(const Matrix<T>& A, const Vector<T>& y, Vector<T>& res)
{
    T sum;
    size_t dim = A.height();
    res = y;

    //last entry
    res[dim-1] = res[dim-1] / A(dim-1, dim-1);

    //remaining entries
    for(int i=dim-2; i>=0; --i)
    {
        //res[i-1] = res[i-1];
        sum = 0;
        for(uint k=i+1; k<dim; ++k)
            sum += A(k,i) * res[k];

        res[i] -= sum;
        res[i] /= A(i,i);
    }
}

//solve system with cholesky decomposition
//requires symmetric posivite definite matrix
//conjugate gradient is more efficient
template <typename T> void NuLLSolve::solveCholesky(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res)
{
    Vector<T> tmp(A.height());
    Matrix<T> chol(A.height());
    NuLLDecomposition::choleskyDecomposition(A, chol);
    forwardElimination(chol, b, tmp);
    backwardSubstitution(chol, tmp, res);
}

//buggy! cross-check!
//solve system with LU-decomposition
//algorithm inspired by gaussian elimination
//applicable to all quadratic matrizes
template <typename T> void NuLLSolve::solveLU(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res)
{
    Vector<T> tmp(A.height());
    Matrix<T> L(A.height());
    Matrix<T> U(A.height());
    NuLLDecomposition::LUDecomposition(A, L, U);
    forwardElimination(L, b, tmp);
    backwardSubstitution(U, tmp, res);
}

//fast method for calculation of equation systems
//needs atmost n = A.dimension() iterations
//slower than thomas algorithm, which works for triagonal systems only
//works only for symmetric positive definite matrizes
template <typename T> void NuLLSolve::conjugateGradient(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res)
{
    size_t dim = A.height();
    Vector<T> rk(dim);
    Vector<T> Apk(dim);
    T alpha, beta, gamma, unscaledAlpha;

    Vector<T> rkN(dim);
    Vector<T> pk(dim);
    Vector<T> x0(dim);

    NuLLTools::MatrixVectorProduct(A, x0, rk);

    for(uint i=0; i<dim; ++i)
        pk[i] = rk[i] = (b[i] - rk[i]);

    for(uint k=0; k<dim; ++k)
    {
        NuLLTools::MatrixVectorProduct(A, pk, Apk);
        gamma = NuLLMath::dotProduct(pk, Apk);

        //check for convergence
        if(!gamma)
            break;

        unscaledAlpha = NuLLMath::dotProduct(rk, rk);
        alpha = unscaledAlpha / gamma;

        for(uint i=0; i<dim; ++i)
        {
            x0[i] = x0[i] + alpha * pk[i];
            rkN[i] = rk[i] - alpha * Apk[i];
        }

        beta = NuLLMath::dotProduct(rkN, rkN) / unscaledAlpha;

        for(uint i=0; i<dim; ++i)
        {
            pk[i] = rkN[i] + beta * pk[i];
            rk[i] = rkN[i];
        }
    }

    res = x0;
}

//cross-check! last position in solution wrong!!
//numerically instable?
//comes very close to solveLU, but is more efficient
//adapted to tridiagonal matrizes; needs less calculations than solveLU
//requires tridiagonal matrix!
template <typename T> void NuLLSolve::ThomasAlgorithm(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res)
{
    size_t dim = A.height();
    Vector<T> m(dim);
    Vector<T> l(dim-1);
    Vector<T> r(dim-1);
    Vector<T> y(dim);
    Vector<T> u(dim);

    //adapted LU decomposition
    m[0] = A(0,0);
    for(uint i=0; i<dim-1; ++i)
    {
        l[i] = A(i,i+1)/m[i];
        m[i+1] = A(i+1,i+1) - l[i] * A(i+1,i);
    }

    //adapted forward elimination
    y[0] = b[0];
    for(uint i=1; i<dim; ++i)
    {
        y[i] = b[i] - l[i-1] * y[i-1];
    }

    //adapted backward substitution
    u[dim-1] = m[dim-1];
    for(int i=dim-2; i>=0; --i)
    {
        u[i] = (y[i] - A(i+1,i) * u[i+1]) / m[i];
    }
}

//calculation of equation systems
//might get trapped in local minimum
template <typename T> void NuLLSolve::gradientDescent(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res, const uint iterations, const T threshold)
{
    uint i = 0;
    size_t dim = A.height();
    Vector<T> uk(dim);
    Vector<T> xk(dim);
    Vector<T> rk(dim);
    T tau, div;

    NuLLTools::MatrixVectorProduct(A, res, rk);
    rk -= b;

    do
    {
        ++i;
        NuLLTools::MatrixVectorProduct(A, rk, uk);
        tau  = NuLLMath::dotProduct(rk, rk);
        div  = NuLLMath::dotProduct(uk, rk);

        if(!div)
            break;

        tau /= div;

        xk *= tau;
        rk *= tau;

        xk -= rk;
        rk -= uk;
    }
    while((tau > threshold) && (i<iterations));

    res = uk;
}

//solves systems for diagonally dominant matrices
//based on matrix splitting
//iterative approach with fast iteration steps
template <typename T> void NuLLSolve::jacobiMethod(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res, const uint iterations, const T threshold)
{
    uint i = 0;
    const size_t dim = A.height();
    T len, lenOld, diff;
    T tmp;

    Vector<T> xn(dim);
    Vector<T> xnp(dim);
    Vector<T> invA(dim);               //diagonal matrix
    xn = xnp = res;

    for(uint i=0; i<dim; ++i)
        invA[i] = 1.0 / A(i,i);

    len = NuLLMath::euclideanNorm(xn);


    do
    {
        ++i;
        lenOld = len;
        xn = xnp;

        for(uint y=0; y<dim; ++y)
        {
            tmp = 0;
            for(uint x=0; x<y; ++x)
                tmp += A(x,y) * xnp[x];

            for(uint x=y+1; x<dim; ++x)
                tmp += A(x,y) * xn[x];

            xnp[y] = invA[y] * (b[y] - tmp);
        }

        len = NuLLMath::euclideanNorm(xnp);
        diff = (len > lenOld)? (len - lenOld) : (lenOld - len);
    }
    while((diff > threshold) && (i<iterations));

    res = xn;
}

//iteratively solves linear systems
//based on matrix splitting
//fast than jacobi method
template <typename T> void NuLLSolve::gaussSeidelMethod(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res, const uint iterations, const T threshold)
{
    uint i = 0;
    const size_t dim = A.height();
    T len, lenOld, diff;
    T tmp;

    Vector<T> xn(dim);
    Vector<T> xnp(dim);
    Vector<T> invA(dim);               //diagonal matrix
    xn = xnp = res;

    for(uint i=0; i<dim; ++i)
        invA[i] = 1.0 / A(i,i);

    len = NuLLMath::euclideanNorm(xn);

    do
    {
        ++i;
        lenOld = len;
        xn = xnp;

        for(uint y=0; y<dim; ++y)
        {
            tmp = 0;
            for(uint x=0; x<y; ++x)
                tmp += A(x,y) * xn[x];

            for(uint x=y+1; x<dim; ++x)
                tmp += A(x,y) * xn[x];

            xnp[y] = invA[y] * (b[y] - tmp);
        }

        len = NuLLMath::euclideanNorm(xnp);
        diff = (len > lenOld)? (len - lenOld) : (lenOld - len);
    }
    while((diff > threshold) && (i<iterations));

    res = xnp;
}

template <typename T> void NuLLSolve::gaussSeidelSOR(const Matrix<T>& A, const Vector<T>&b, Vector<T>& res, const T omega, const uint iterations = 100, const T threshold = 0.00001)
{
    uint i = 0;
    const size_t dim = A.height();
    const T omegaM = 1. - omega;
    T len, lenOld, diff;
    T tmp;

    Vector<T> xn(dim);
    Vector<T> xnp(dim);
    Vector<T> invA(dim);               //diagonal matrix
    xn = xnp = res;

    for(uint i=0; i<dim; ++i)
        invA[i] = omega / A(i,i);

    len = NuLLMath::euclideanNorm(xn);

    do
    {
        ++i;
        lenOld = len;
        xn = xnp;

        for(uint y=0; y<dim; ++y)
        {
            tmp = 0;
            for(uint x=0; x<y; ++x)
                tmp += A(x,y) * xn[x];

            for(uint x=y+1; x<dim; ++x)
                tmp += A(x,y) * xn[x];

            xnp[y] = invA[y] * (b[y] - tmp);
            xnp[y] += omegaM * xn[y];
        }

        len = NuLLMath::euclideanNorm(xnp);
        diff = (len > lenOld)? (len - lenOld) : (lenOld - len);
    }
    while((diff > threshold) && (i<iterations));

    res = xnp;
}

#endif
