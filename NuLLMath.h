#ifndef NULLMATH_H
#define NULLMATH_H

#include "Utils.h"
#include "NuLLDecomposition.inl"
#include "NuLLSolve.inl"
#include "Matrix.h"
#include "Vector.h"

/*
 * various mathematical methods for vectors and matrices
 * used by other NuLL headers
 * requires NuLLDecomposition.h
 *
 */

namespace NuLLMath
{
    const double pi = 2.0 * asinf(1.0);
    const double e  = 2.7182818284590452353602874713526624977572470936999595;

    template <typename T> void transpose(const Matrix<T>& mtx, Matrix<T>& dst);

    template <typename T> T determinant(const Matrix<T>& mtx);
    template <typename T> T determinant2x2(const Matrix<T>& mtx);
    template <typename T> T determinant3x3(const Matrix<T>& mtx);

    template <typename T> void invert(const Matrix<T>& mtx, Matrix<T>& dst);
    template <typename T> void invert2x2(const Matrix<T>& mtx, Matrix<T>& dst);
    template <typename T> void invert3x3(const Matrix<T>& mtx, Matrix<T>& dst);

    template <typename T> bool diagonallyDominant(const Matrix<T>& mtx);

    template <typename T> void dyadicProduct(const Vector<T>& a, const Vector<T>& b, Matrix<T>& dst);
    template <typename T> void dyadicProduct(const Vector<T>& a, Matrix<T>& dst);
    template <typename T> void kroneckerProduct(const Matrix<T>& A, Matrix<T>& B, Matrix<T>& dst);

    template <typename T> T dotProduct(const Vector<T>& x, const Vector<T>& y);
    template <typename T> T euclideanNorm(const Vector<T>& x);

    template <typename T> T columnSumNorm(const Matrix<T>& mtx);
    template <typename T> T rowSumNorm(const Matrix<T>& mtx);
    template <typename T> T frobeniusNorm(const Matrix<T>& mtx);
}

#endif // NULLMATH_H
