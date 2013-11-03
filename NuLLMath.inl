#ifndef NULLMATH_INL
#define NULLMATH_INL

#include "NuLLMath.h"

//transpose matrix
//useless for symmatrix (doh!)
template <typename T> void NuLLMath::transpose(const Matrix<T>& mtx, Matrix<T>& dst)
{
    size_t width = mtx.width();
    size_t height = mtx.height();

    for(uint y=0; y<height; ++y)
        for(uint x=0; x<width; ++x)
            dst(x,y) = mtx(y,x);
}

//hardcoded algorithm for calculation of determinant (2x2)
template <typename T> T NuLLMath::determinant2x2(const Matrix<T>& mtx)
{
    return (mtx(0,0) * mtx(1,1) - mtx(0,1) * mtx(1,0));
}

//hardcoded algorithm for calculation of determinant (3x3)
template <typename T> T NuLLMath::determinant3x3(const Matrix<T>& mtx)
{
    T res = 0;
    res += mtx(0,0) * mtx(1,1) * mtx(2,2);
    res += mtx(0,1) * mtx(1,2) * mtx(2,0);
    res += mtx(0,2) * mtx(1,0) * mtx(2,1);
    res -= mtx(0,2) * mtx(1,1) * mtx(2,0);
    res -= mtx(0,1) * mtx(1,0) * mtx(2,2);
    res -= mtx(0,0) * mtx(1,2) * mtx(2,1);
    return res;
}

//TODO: check!
//simple algorithm for determinant calculation
//uses LU decomposition for big matrizes
//switches to special cases for mtx.dimension() < 4
template <typename T> T NuLLMath::determinant(const Matrix<T>& mtx)
{
    size_t dim = mtx.height();

    if(dim == 1)
        return mtx(0,0);

    if(dim == 2)
        return (determinant2x2(mtx));

    if(dim == 3)
        return (determinant3x3(mtx));

    //calculation with LU decomposition
    //det(A) = det(L*U) = det(L) * det(U) = det(U) = prod(U(i,i))
    //multiply entries of U matrix only since diagonals of L matrix are one
    T res = 1;
    PowerMatrix<T> L(dim);
    PowerMatrix<T> U(dim);
    NuLLDecomposition::LUDecomposition(mtx, L, U);

    for(uint i=0; i<dim; ++i)
        res *= U(i,i);

    return res;
}

template<typename T> void NuLLMath::invert2x2(const Matrix<T>& mtx, Matrix<T>& dst)
{
    T factor = determinant2x2(mtx);
    dst(0,0) = mtx(1,1);
    dst(0,1) = -mtx(0,1);
    dst(1,0) = -mtx(1,0);
    dst(1,1) = mtx(0,0);
    dst /= factor;
}
//TODO: implement
template<typename T> void NuLLMath::invert3x3(const Matrix<T>& mtx, Matrix<T>& dst)
{
    T factor = determinant3x3(mtx);
    dst(0,0) = mtx(1,1) * mtx(2,2) - mtx(1,2) * mtx(2,1);
    dst(0,1) = mtx(2,1) * mtx(0,2) - mtx(0,1) * mtx(2,2);
    dst(0,2) = mtx(0,1) * mtx(1,2) - mtx(1,1) * mtx(0,2);
    dst(1,0) = mtx(3,0) * mtx(1,2) - mtx(0,1) * mtx(2,2);
    dst(1,1) = mtx(0,0) * mtx(2,2) - mtx(0,2) * mtx(2,0);
    dst(1,2) = mtx(2,0) * mtx(0,1) - mtx(0,0) * mtx(2,1);
    dst(2,0) = mtx(0,1) * mtx(2,1) - mtx(1,1) * mtx(0,2);
    dst(2,1) = mtx(1,0) * mtx(0,2) - mtx(0,0) * mtx(1,2);
    dst(2,2) = mtx(0,0) * mtx(1,1) - mtx(0,1) * mtx(1,0);
    dst /= factor;
}

//TODO: implement and check!
//uses QR decomposition for fast inversion:
// A^-1 = (QR)^-1 = (R^-1)(Q^-1) = (R^-1)(transpose(Q))
//hard-coded cases for 2x2 and 3x3 matrizes
template<typename T> void NuLLMath::invert(const Matrix<T>& mtx, Matrix<T>& dst)
{
    size_t dim = mtx.width();

    if(dim == 1)
    {
        dst(0,0) = (1/mtx(0,0));
        return;
    }
    if(dim == 2)
    {
        invert2x2(mtx, dst);
        return;
    }
    if(dim == 3)
    {
        invert3x3(mtx, dst);
        return;
    }

    //LU inversion
    PowerMatrix<T> L(mtx.width(), mtx.height());
    PowerMatrix<T> U(mtx.width(), mtx.height());
    //NuLLDecomposition::LUDecomposition(mtx, L, U);
}

//quick check if matrix is positive definite
//uses criteria for strict diagonal dominant matrizes
template <typename T> bool NuLLMath::diagonallyDominant(const Matrix<T>& mtx)
{
    T colSum;
    size_t height = mtx.height();
    size_t width = mtx.width();

    for(uint x=0; x<width; ++x)
    {
        colSum = 0;
        for(uint y=0; y<height; ++y)
        {
            if(x==y)
                continue;

            colSum += abs(mtx(x,y));
        }
        std::cout << colSum << "/" << mtx(x,x) << std::endl;
        if(colSum >= abs(mtx(x,x)))
            return false;
    }
    return true;
}

//dyadic product of two vectors
template <typename T> void NuLLMath::dyadicProduct(const Vector<T>& a, const Vector<T>& b, Matrix<T>& dst)
{
    for(uint x=0; x<a.size(); ++x)
        for(uint y=0; y<b.size(); ++y)
            dst(x,y) = a[x] * b[y];
}

//dyadic product of a vector with itself
//faster than calling the naive function via dyadicProduct(a, a)
template <typename T> void NuLLMath::dyadicProduct(const Vector<T>& a, Matrix<T>& dst)
{
    for(uint x=0; x<a.size(); ++x)
        for(uint y=x; y<a.size(); ++y)
            dst(y,x) = dst(x,y) = a[x] * a[y];
}

//kronecker product
template <typename T> void NuLLMath::kroneckerProduct(const Matrix<T>& A, Matrix<T>& B, Matrix<T>& dst)
{
    for(uint ya=0; ya<A.height(); ++ya)
    {
        for(uint xa=0; xa<A.width(); ++xa)
        {
            for(uint yb=0; yb<B.height(); ++yb)
            {
                for(uint xb=0; xb<B.width(); ++xb)
                {
                    int x = xb+xa*B.width();
                    int y = yb+ya*B.height();
                    //std::cout << x << "/" << y << std::endl;
                    dst(x,y) = A(xa,ya) * B(xb,yb);
                }
            }
        }
    }
}

//simple dot product
//be aware that vectors need to have same dimensions
//remember: returns zero if vectors are orthogonal
template <typename T> T NuLLMath::dotProduct(const Vector<T>& x, const Vector<T>& y)
{
    T res = 0;
    if(x.size() != y.size())
        return res;

    for(uint i=0; i<x.size(); ++i)
        res += x[i] * y[i];

    return res;
}

//euclidian norm
template <typename T> T NuLLMath::euclideanNorm(const Vector<T>& x)
{
    return sqrt(dotProduct(x, x));
}

//column sum norm
template <typename T> T NuLLMath::columnSumNorm(const Matrix<T>& mtx)
{
    size_t width = mtx.width();
    size_t height = mtx.height();

    T sum;
    T res = 0;
    for(uint x=0; x<width; ++x)
    {
        sum = 0;
        for(uint y=0; y<height; ++y)
        {
            sum += mtx(x,y);
        }

        if(sum > res)
            res = sum;
    }
    return res;
}

//row sum norm
template <typename T> T NuLLMath::rowSumNorm(const Matrix<T>& mtx)
{
    size_t width = mtx.width();
    size_t height = mtx.height();

    T sum;
    T res = 0;
    for(uint x=0; x<width; ++x)
    {
        sum = 0;
        for(uint y=0; y<height; ++y)
        {
            sum += mtx(y,x);
        }

        if(sum > res)
            res = sum;
    }
    return res;
}

//frobenius norm
//equivalent to euclidean norm of vectors
template <typename T> T NuLLMath::frobeniusNorm(const Matrix<T>& mtx)
{
    size_t width = mtx.width();
    size_t height = mtx.height();

    T res = 0;
    for(uint x=0; x<width; ++x)
    {
        for(uint y=0; y<height; ++y)
        {
            res += mtx(x,y) * mtx(x,y);
        }
    }
    return sqrt(res);
}

#endif
