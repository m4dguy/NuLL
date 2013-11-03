#ifndef NULLEIGEN_H
#define NULLEIGEN_H

#include "Vector.h"
#include "RectMatrix.h"
#include "NuLLDecomposition.h"

/*
 * methods for eigenvalue and eigenvector calculation
 *
 */

namespace NuLLEigen
{
    template<typename T> inline void powerIteration(const Matrix<T>& mtx, Vector<T>& dstVec, T& dstVal, T threshold = 0.00001);
    template<typename T> inline void QRAlgorithm(const Matrix<T>& A, RectMatrix<T>& dstVecs, Vector<T>& dstVals, T threshold = 0.00001);
}

#endif // NULLEIGEN_H
