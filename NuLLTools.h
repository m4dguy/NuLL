#ifndef NULLTOOLS_H
#define NULLTOOLS_H

#include "Utils.h"
#include "Vector.h"
#include "Matrix.h"

namespace NuLLTools
{
    template <typename T> void copyMatrix(const Matrix<T>& src, Matrix<T>& dst);
    template <typename T> void getSegment(const Matrix<T>& src, Matrix<T>& dst, int startX=0, int startY=0, int endX=1, int endY=1);
	template <typename T> void pasteAt(Matrix<T>& mtx, const Matrix<T>& other, int dstX=0, int dstY=0);

    template <typename T> void getColumn(const Matrix<T>& src, Vector<T>& dst, uint col);
    template <typename T> void getRow(const Matrix<T>& src, Vector<T>& dst, uint row);

	template <typename T> T maxValue(const Matrix<T>& mtx);
	template <typename T> T minValue(const Matrix<T>& mtx);
	template <typename T> T average(const Matrix<T>& mtx);
}

#endif // NULLTOOLS_H
