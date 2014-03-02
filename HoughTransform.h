#ifndef HOUGHTRANSFORM_H
#define HOUGHTRANSFORM_H

#include <vector>

#include "Matrix.h"

namespace HoughTransform
{
	typedef struct
	{
		uint x;
		uint y;
		uint radius;
		double score;
	}HoughCircle;

	template<typename T> void houghCircles(const Matrix<T>& img, std::vector<HoughCircle>& circles, double threshold=100, int minRad=1, int maxRad=50);
	template<typename T> void plotCircles(const std::vector<HoughCircle>& circles, Matrix<T>& img, uint radius=1);

	void mergeCircles(std::vector<HoughCircle>& circles, uint maxDist=10);
	template<typename T> void histogramInCircle(const Matrix<T>& mtx, const HoughCircle& circ, std::vector<uint>& hist);
}

#endif