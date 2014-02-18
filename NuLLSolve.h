#ifndef NULLSOLVE_H
#define NULLSOLVE_H

#include "Vector.h"
#include "NuLLDecomposition.inl"
#include "NuLLTools.inl"

/*
 * different methods for solving linear systems of equations
 * some methods require matrix decompositions
 *
 */

namespace NuLLSolve
{
    template<typename T> void forwardElimination(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);
    template<typename T> void backwardSubstitution(const Matrix<T>& A, const Vector<T>& y, Vector<T>& res);

    template <typename T> void solveCholesky(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);
    template <typename T> void solveLU(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);

    template <typename T> void conjugateGradient(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);
    template <typename T> void ThomasAlgorithm(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);

    template <typename T> void jacobiMethod(const Matrix<T> A, const Vector<T>& b, Vector<T>& res, T threshold = 0.00001);

}

#endif // NULLSOLVE_H
