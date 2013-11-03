#ifndef NULLTOOLS_INL
#define NULLTOOLS_INL

#include "NuLLTools.h"

//can be used for matrix conversion
//copy content of matrix src to matrix dst
template <typename T> void NuLLTools::copyMatrix(const Matrix<T>& src, Matrix<T>& dst)
{
    size_t width = src.width();
    size_t height = src.height();
    for(uint y=0; y<height; ++y)
        for(uint x=0; x<width; ++x)
            dst(x,y) = src(x,y);
}

//copy segment of matrix to dst
//uses absolute coordinates
template <typename T> void NuLLTools::getSegment(const Matrix<T>& src, Matrix<T>& dst, int startX, int startY, int endX, int endY)
{

    if(endX < startX)
        std::swap(startX, endX);

    if(endY < startY)
        std::swap(startY, endY);

    for(int y=startY; y<endY; ++y)
        for(int x=startX; x<endX; ++x)
            dst(x,y) = src.getMirrored(x,y);
}

//pastes matrix "other" into matrix "mtx"
template <typename T> void NuLLTools::pasteAt(Matrix<T>& mtx, const Matrix<T>& other, int dstX, int dstY)
{
	size_t width = other.width();
	size_t height = other.height();
	for(uint y=0; y<width; ++y)
		for(uint x=0; x<height; ++x)
			mtx(x+dstX, y+dstY) = other(x,y);
}

//copy matrix column to vector
template <typename T> void NuLLTools::getColumn(const Matrix<T>& src, Vector<T>& dst, uint col)
{
    size_t dim = src.height();
    for(uint y=0; y<dim; ++y)
        dst[y] = src(col, y);

}

//copy matrix row to vector
template <typename T> void NuLLTools::getRow(const Matrix<T>& src, Vector<T>& dst, uint row)
{
    size_t dim = src.width();
    for(uint x=0; x<dim; ++x)
        dst[x] = src(x, row);

}
template <typename T> T NuLLTools::maxValue(const Matrix<T>& mtx)
{
	T maxV = mtx(0,0);
	size_t width = mtx.width();
	size_t height = mtx.height();

	for(uint y=0; y<height; ++y)
		for(uint x=0; x<width; ++x)
			maxV = max(maxV, mtx(x,y));

	return maxV;
}

template <typename T> T NuLLTools::minValue(const Matrix<T>& mtx)
{
	T minV = mtx(0,0);
	size_t width = mtx.width();
	size_t height = mtx.height();

	for(uint y=0; y<height; ++y)
		for(uint x=0; x<width; ++x)
			minV = min(minV, mtx(x,y));

	return minV;
}

template <typename T> T NuLLTools::average(const Matrix<T>& mtx)
{
	T avg = 0;
	size_t width = mtx.width();
	size_t height = mtx.height();

	for(uint y=0; y<height; ++y)
		for(uint x=0; x<width; ++x)
			avg += mtx(x,y);

	avg /= width * height;
	return avg;
}

#endif
