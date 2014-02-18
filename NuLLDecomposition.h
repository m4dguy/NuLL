#ifndef NULLDECOMPOSITION_H
#define NULLDECOMPOSITION_H

#include <math.h>

#include "Utils.h"
#include "NuLLMath.inl"
#include "NuLLTools.inl"

#include "Vector.h"
#include "Matrix.h"

/*
 * different matrix decomposition methods
 *
 */

namespace NuLLDecomposition
{
    template <typename T> void choleskyDecomposition(const Matrix<T>& mtx, Matrix<T>& dst);
    template <typename T> void LUDecomposition(const Matrix<T>& mtx, Matrix<T>& dstL, Matrix<T>& dstU);
    template <typename T> void QRDecomposition(const Matrix<T>& mtx, Matrix<T>& dstQ, Matrix<T>& dstR);

}

#endif // NULLDECOMPOSITION_H
