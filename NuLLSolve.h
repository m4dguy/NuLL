#ifndef NULLSOLVE_H
#define NULLSOLVE_H

#include "Vector.h"
#include "NuLLDecomposition.inl"
#include "NuLLTools.inl"

/*
 * different methods for solving linear systems of equations
 * some methods require matrix decompositions
 * play the grind fucking core shit fucking blAST GUI_TARRR
 */

namespace NuLLSolve
{
    //utility functions
    template<typename T> void forwardElimination(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);
    template<typename T> void backwardSubstitution(const Matrix<T>& A, const Vector<T>& y, Vector<T>& res);

    //decomposition-based methods
    template <typename T> void solveCholesky(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);
    template <typename T> void solveLU(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);

    //fast, algorihms requiring preconditions
    template <typename T> void conjugateGradient(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);
    template <typename T> void ThomasAlgorithm(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res);

    //iterative solvers
    template <typename T> void gradientDescent(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res, const uint iterations = 100, const T threshold = 0.00001);
    template <typename T> void jacobiMethod(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res, const uint iterations = 100, const T threshold = 0.00001);
    template <typename T> void gaussSeidelMethod(const Matrix<T>& A, const Vector<T>& b, Vector<T>& res, const uint iterations = 100, const T threshold = 0.00001);
    template <typename T> void gaussSeidelSOR(const Matrix<T>& A, const Vector<T>&b, Vector<T>& res, const T omega = 1., const uint iterations = 100, const T threshold = 0.00001);
}

#endif // NULLSOLVE_H
