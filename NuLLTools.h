#ifndef NULLTOOLS_H
#define NULLTOOLS_H

#include <math.h>

#include "Utils.h"
#include "Vector.h"
#include "Matrix.h"

namespace NuLLTools
{
	//copying of entries
	template <typename T, typename U> void copyMatrix(const Matrix<T>& src, Matrix<U>& dst);
	template <typename T> void getSegment(const Matrix<T>& src, Matrix<T>& dst, int startX=0, int startY=0, int endX=1, int endY=1);
	template <typename T> void pasteAt(const Matrix<T>& mtx, Matrix<T>& dst, int dstX=0, int dstY=0);

	//getting rows and columns
	template <typename T> void getColumn(const Matrix<T>& src, Vector<T>& dst, uint col);
	template <typename T> void getRow(const Matrix<T>& src, Vector<T>& dst, uint row);

	//elementwise operations
	template <typename T> void elementwiseMultiplication(const Matrix<T>& mtx1, const Matrix<T> mtx2, Matrix<T>& dst);
	template <typename T> void elementwiseAddition(const Matrix<T>& mtx1, const Matrix<T> mtx2, Matrix<T>& dst);
	template <typename T> void elementwiseEquals(const Matrix<T>& mtx1, const Matrix<T> mtx2, Matrix<T>& dst);

	//utility
	template <typename T> T maxValue(const Matrix<T>& mtx);
	template <typename T> T minValue(const Matrix<T>& mtx);
	template <typename T> T average(const Matrix<T>& mtx);
	template <typename T> void abs(const Matrix<T>& mtx, Matrix<T>& dst);
	template <typename T> void sqrt(const Matrix<T>& mtx, Matrix<T>& dst);
}

#endif // NULLTOOLS_H
