#ifndef NULLPROCESSING_H
#define NULLPROCESSING_H

#include <algorithm>
#include <vector>
#include <math.h>

#include "Utils.h"
#include "NuLLTools.inl"

/*
 *
 * methods for data and signal processing
 * includes:    point operations
 *              morphological filters
 *              derivative filters (first order, second order)
 *              kernel generators
 *              convolution
 *
 */

namespace NuLLProcessing
{
	//kernel operations
    template <typename T> void normalize(Matrix<T>& mtx);
    template <typename T> void convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Matrix<T>& kernel);

	//kernel generators
	template <typename T> void identityKernel(Matrix<T>& dst, int radius = 1);
	template <typename T> void pillboxKernel(Matrix<T>& dst, int radius = 1);
    template <typename T> void gaussianKernel(Matrix<T>& dst, int radius = 1, double sigma = 1.0);

	//simple blur techniques
    template <typename T> void pillboxBlur(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
    template <typename T> void gaussianBlur(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1, double sigma = 1.0);

	//derivative filters
    template <typename T> void firstDerivative(const Matrix<T>& mtx, Matrix<T>& dst);
    template <typename T> void secondDerivative(const Matrix<T>& mtx, Matrix<T>& dst);

	//others
	template <typename T> void localVariance(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void localCurvature(const Matrix<T>& mty, Matrix<T>& dst, int radius = 1);

	//morphological filters
    template <typename T> void medianFilter(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1, float percentile = .5f);
    template <typename T> void dilation(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
    template <typename T> void erosion(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
    template <typename T> void opening(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
    template <typename T> void closing(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void whiteTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void blackTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);
	template <typename T> void selfdualTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius = 1);

	//point operations
	template <typename T> void gammaCorrection(const Matrix<T>& mtx, Matrix<T>& dst, double gamma = 1.2);
	template <typename T> void snap(const Matrix<T>& mtx, Matrix<T>& dst, T min = 0, T max = 255);
	template <typename T> void threshold(const Matrix<T>& mtx, Matrix<T>& dst, double threshold, double gmin = 0, double gmax = 255);
	template <typename T> void doubleThreshold(const Matrix<T>& mtx, Matrix<T>& dst, double thresholdLower, double thresholdUpper, double gmin = 0, double gmax = 255);
	template <typename T> void automatedThreshold(const Matrix<T>& mtx, Matrix<T>& dst, double gmin = 0, double gmax = 255);
	template <typename T> void logDynamicCompression(const Matrix<T>& mtx, Matrix<T>& dst, double c = 0);
	template <typename T> void affineRescale(const Matrix<T>& mtx, Matrix<T>& dst, double minVal = 0, double maxVal = 255);
	template <typename T> void affineTransform(const Matrix<T>& mtx, Matrix<T>& dst, double a = 1, double b = 0);

	//noise mask generators
	//use for additive or multiplicative noise
	//use of "snap" after addition/ multiplication recommended to conserve range [0,255]
	template <typename T> void uniformNoiseMask(Matrix<T>& dst, T a, T b, double percentage = 0.5);
	template <typename T> void gaussianNoiseMask(Matrix<T>& dst, double sigma, double my, double percentage = 0.5);
	template <typename T> void impulseNoiseMask(Matrix<T>& dst, double percentage = 0.5);
}

#endif // NULLPROCESSING_H
