#ifndef NULLTOOLS_H
#define NULLTOOLS_H

#include <math.h>
#include <stdarg.h>

#include "Utils.h"
#include "Vector.h"
#include "Matrix.h"

namespace NuLLTools
{
	//copying of entries
	template <typename T, typename U> void copyMatrix(const Matrix<T>& src, Matrix<U>& dst);
	template <typename T> void getSegment(const Matrix<T>& src, Matrix<T>& dst, const int startX=0, const int startY=0, const int endX=1, const int endY=1);
	template <typename T> void pasteAt(const Matrix<T>& mtx, Matrix<T>& dst, const int dstX=0, const int dstY=0);
    template <typename T> void makeVector(Vector<T>& dst, int elements, ...);
    template <typename T> void makeMatrix(Matrix<T>& dst, int elements, ...);

	//getting rows and columns
	template <typename T> void getColumn(const Matrix<T>& src, Vector<T>& dst, const uint col);
	template <typename T> void getRow(const Matrix<T>& src, Vector<T>& dst, const uint row);

    //Multiplications
    template <typename T> void MatrixVectorProduct(const Matrix<T>& mtx, const Vector<T>& vec, Vector<T>& dst);
    template <typename T> void MatrixMatrixProduct(const Matrix<T>& mtx1, const Matrix<T>& mtx2, Matrix<T>& dst);

	//elementwise operations
	template <typename T> void elementwiseMultiplication(const Matrix<T>& mtx1, const Matrix<T> mtx2, Matrix<T>& dst);
	template <typename T> void elementwiseAddition(const Matrix<T>& mtx1, const Matrix<T> mtx2, Matrix<T>& dst);
	template <typename T> void elementwiseEquals(const Matrix<T>& mtx1, const Matrix<T> mtx2, Matrix<T>& dst);

	//statistics
	template <typename T> std::pair<T,T> minMaxValue(const Matrix<T>& mtx);
	template <typename T> T maxValue(const Matrix<T>& mtx);
	template <typename T> T minValue(const Matrix<T>& mtx);
	template <typename T> T mean(const Matrix<T>& mtx);
	template <typename T> T variance(const Matrix<T>& mtx);

	//utility
	template <typename T> void mAbs(const Matrix<T>& mtx, Matrix<T>& dst);
	template <typename T> void mSqrt(const Matrix<T>& mtx, Matrix<T>& dst);

}

#endif // NULLTOOLS_H
