#ifndef NULLPROCESSING_INL
#define NULLPROCESSING_INL

#include "NuLLProcessing.h"

//discrete convolution of matrix (2D convolution)
template<typename T> void NuLLProcessing::convolve(const Matrix<T>& mtx, Matrix<T>& dst, const Matrix<T>& kernel)
{
    T val;
    size_t width = mtx.width();
    size_t height = mtx.height();

    size_t kWidth = kernel.width();
    size_t kHeight = kernel.height();

    size_t offsetX = kernel.width() / 2;
    size_t offsetY = kernel.height() / 2;

    for(uint y=0; y<height; ++y)
    {
        for(uint x=0; x<width; ++x)
        {
            val = 0.0f;
            for(uint ky=0; ky<kHeight; ++ky)
            {
                for(uint kx=0; kx<kWidth; ++kx)
                {
                    val += kernel(kx, ky) * mtx.getMirrored(x+(offsetX-kx), y+(offsetY-ky));
                }
            }
            dst(x,y) = val;
        }
    }
}

//normalizes matrix to make sum of entries equal 1
//useful for kernel normalization
template<typename T> void NuLLProcessing::normalize(Matrix<T>& mtx)
{
    T factor = 0;
    size_t height = mtx.height();
    size_t width = mtx.width();

    for(uint y=0; y<height; ++y)
        for(uint x=0; x<width; ++x)
            factor += mtx(x,y);

    mtx /= factor;
}

//creates identity kernel for convolution
//useful for kernel design
template <typename T> void NuLLProcessing::identityKernel(Matrix<T>& dst, int radius)
{
	const uint dim = (2 * radius) + 1;
	for(uint y=0; y<dim; ++y)
		for(uint x=0; x<dim; ++x)
			dst(x,y) = 0;

	dst(radius, radius) = 1.0;
}

template <typename T> void NuLLProcessing::pillboxKernel(Matrix<T>& dst, int radius)
{
    const int dim = (2 * radius) + 1;
    int radSq = radius * radius;
    T scale = 0;

    //kernel calculation
    for(int y=0; y<dim; ++y)
    {
        for(int x=0; x<dim; ++x)
        {
            if(((radius-x)*(radius-x) + (radius-y) * (radius-y)) <= radSq)
            {
                dst(x,y) = 1.0;
                ++scale;
            }
        }
    }
    dst /= scale;
}

//creates a gaussian kernel of given size
//sigma/ variance parameter only has notable effect for large kernels
template <typename T> void NuLLProcessing::gaussianKernel(Matrix<T>& dst, int radius, double sigma)
{
    const uint dim = (2 * radius) + 1;
	const double pi = 2.0 * asinf(1.0);

    T div = 2 * pi * sigma * sigma;
    T kx, ky, res, scale;
    scale = 0;

    //kernel calculation
    for(uint y=0; y<dim; ++y)
    {
        for(uint x=0; x<dim; ++x)
        {
            kx = (x - radius) * (x - radius);
            ky = (y - radius) * (y - radius);
            res = -(ky + kx);
            res = exp(res / (2 * sigma * sigma));
            res /= div;
            dst(x,y) = res;

            scale += res;
        }
    }

	//normalization: sum of entries only roughly == 1
    dst /= scale;
}

//smoothing with pillbox kernel
template <typename T> void NuLLProcessing::pillboxBlur(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    Matrix<T> kernel((2 * radius) + 1);
    pillboxKernel(kernel, radius);
    convolve(mtx, dst, kernel);
}

//gaussian smoothing for matrices
//update me: use separation theorem for fast gaussian convolution!
template <typename T> void NuLLProcessing::gaussianBlur(const Matrix<T>& mtx, Matrix<T>& dst, int radius, double sigma)
{
    Matrix<T> kernel(radius+radius+1);
    gaussianKernel(kernel, radius, sigma);
    convolve(mtx, dst, kernel);
}

//first order derivative (central differential quotient) for matrices
template <typename T> void NuLLProcessing::firstDerivative(const Matrix<T>& mtx, Matrix<T>& dst)
{
    int width = mtx.width();
    int height = mtx.height();
    T dx, dy, gradMag;

    for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            dx = mtx.getMirrored(x+1,y);
            dx -= mtx.getMirrored(x-1, y);
            dx /= 2;

            dy = mtx.getMirrored(x, y-1);
            dy -= mtx.getMirrored(x, y+1);
            dy /= 2;

            gradMag = sqrt(dx*dx+dy*dy);
            dst(x,y) = gradMag;
        }
    }
}

//second order derivative (central differential quotient) for matrices
template <typename T> void NuLLProcessing::secondDerivative(const Matrix<T>& mtx, Matrix<T>& dst)
{
    int width = mtx.width();
    int height = mtx.height();
    T dx, dy, gradMag;

    for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            dx = mtx.getMirrored(x-1, y);
            dx -= 2 * mtx.getMirrored(x, y);
            dx += mtx.getMirrored(x+1, y);

            dy = mtx.getMirrored(x, y+1);
            dy -= 2 * mtx.getMirrored(x, y);
            dy += mtx.getMirrored(x, y-1);

            gradMag = sqrt(dx*dx+dy*dy);
            dst(x,y) = gradMag;
        }
    }
}

//result similar to derivative filter
//however, bigger neighborhoods can be respected
template <typename T> void NuLLProcessing::localVariance(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	const int radsq = radius * radius;
	int width = mtx.width();
    int height = mtx.height();
	T mean, var;

	std::vector<T> neighbors;
	neighbors.reserve((radius+1)*(radius+1));

	for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            for(int ky=-radius; ky<=radius; ++ky)
            {
                for(int kx=-radius; kx<=radius; ++kx)
                {
                    //disc-shaped neighborhood
                    if(kx*kx + ky*ky > radsq)
                        continue;

                    neighbors.push_back(mtx.getMirrored(x+kx, y+ky));
                }
            }

			mean = 0;
			for(uint i=0; i<neighbors.size(); ++i)
				mean += neighbors[i];

			mean /= neighbors.size();

			var = 0;
			for(uint i=0; i<neighbors.size(); ++i)
				var += (neighbors[i] - mean) * (neighbors[i] - mean);

			var /= neighbors.size();
			dst(x,y) = sqrt(var);

            neighbors.clear();
        }
    }
}

//TODO: implement me!
//curvature filter
template <typename T> void NuLLProcessing::localCurvature(const Matrix<T>& mty, Matrix<T>& dst, int radius)
{


}

//median filtering for matrices
template <typename T> void NuLLProcessing::medianFilter(const Matrix<T>& mtx, Matrix<T>& dst, int radius, float percentile)
{
    const int radsq = radius * radius;
    int width = mtx.width();
    int height = mtx.height();

    std::vector<T> neighbors;
	neighbors.reserve((radius+1)*(radius+1));

    for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            for(int ky=-radius; ky<=radius; ++ky)
            {
                for(int kx=-radius; kx<=radius; ++kx)
                {
                    //disc-shaped neighborhood
                    if(kx*kx + ky*ky > radsq)
                        continue;

                    neighbors.push_back(mtx.getMirrored(x+kx, y+ky));
                }
            }
            std::sort(neighbors.begin(), neighbors.end());
            dst(x,y) = neighbors[(uint)((neighbors.size()-1) * percentile)];
            neighbors.clear();
        }
    }
}

//dilation for matrices
template <typename T> void NuLLProcessing::dilation(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    int width = mtx.width();
    int height = mtx.height();

    T sup;
    std::vector<T> neighbors;
    neighbors.reserve((radius+1)*(radius+1));
    const int radsq = radius * radius;

    for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            for(int ky=-radius; ky<=radius; ++ky)
            {
                for(int kx=-radius; kx<=radius; ++kx)
                {
                    //disc-shaped neighborhood
                    if(kx*kx + ky*ky > radsq)
                        continue;

                    neighbors.push_back(mtx.getMirrored(x+kx, y+ky));
                }
            }

            sup = neighbors[0];
            for(uint i=1; i<neighbors.size(); ++i)
            {
                if(neighbors[i] < sup)
                    sup = neighbors[i];
            }
            dst(x,y) = sup;
            neighbors.clear();
        }
    }
}

//erosion for matrices
template <typename T> void NuLLProcessing::erosion(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    int width = mtx.width();
    int height = mtx.height();

    T inf;
    std::vector<T> neighbors;
    neighbors.reserve((radius+1)*(radius+1));
    const int radsq = radius * radius;

    for(int y=0; y<height; ++y)
    {
        for(int x=0; x<width; ++x)
        {
            for(int ky=-radius; ky<=radius; ++ky)
            {
                for(int kx=-radius; kx<=radius; ++kx)
                {
                    //disc-shaped neighborhood
                    if(kx*kx + ky*ky > radsq)
                        continue;

                    neighbors.push_back(mtx.getMirrored(x+kx, y+ky));
                }
            }

            inf = neighbors[0];
            for(uint i=1; i<neighbors.size(); ++i)
            {
                if(neighbors[i] > inf)
                    inf = neighbors[i];
            }
            dst(x,y) = inf;
            neighbors.clear();
        }
    }
}

//opening for matrices
template <typename T> void NuLLProcessing::opening(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    Matrix<T> tmp(mtx.width(), mtx.height());
    erosion(mtx, tmp, radius);
    dilation(tmp, dst, radius);
}

//closing for matrices
template <typename T> void NuLLProcessing::closing(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
    Matrix<T> tmp(mtx.width(), mtx.height());
    dilation(mtx, tmp, radius);
    erosion(tmp, dst, radius);
}

//white tophat operation for matrices
template <typename T> void NuLLProcessing::whiteTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
    opening(mtx, tmp, radius);
    dst = mtx;
	dst -= tmp;
}

//black tophat operation for matrices
template <typename T> void NuLLProcessing::blackTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
	closing(mtx, dst, radius);
	dst -= mtx;
}

template <typename T> void NuLLProcessing::selfdualTopHat(const Matrix<T>& mtx, Matrix<T>& dst, int radius)
{
	Matrix<T> tmp(mtx.width(), mtx.height());
	blackTopHat(mtx, dst, radius);
	whiteTopHat(mtx, tmp, radius);
	dst += tmp;
}

//fixes values which are below/ over given range
template <typename T> void NuLLProcessing::snap(const Matrix<T>& mtx, Matrix<T>& dst, T min, T max)
{
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = std::min<T>(std::max<T>(mtx(x,y), min), max);
}

//filters out entries below threshold
template <typename T> void NuLLProcessing::threshold(const Matrix<T>& mtx, Matrix<T>& dst, double threshold, double gmin, double gmax)
{
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = (mtx(x,y) >= threshold)? gmax : gmin;
}

//filter entries outside range
template <typename T> void NuLLProcessing::doubleThreshold(const Matrix<T>& mtx, Matrix<T>& dst, double thresholdLower, double thresholdUpper, double gmin, double gmax)
{
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = (mtx(x,y) >= thresholdLower && mtx(x,y) <= thresholdUpper)? gmax : gmin;
}

//automatically calculate threshold value and filter out entries
//uses otsu's method
template <typename T> void NuLLProcessing::automatedThreshold(const Matrix<T>& mtx, Matrix<T>& dst, double gmin, double gmax)
{
	double threshold = 0;

	double curVar = 0;
	double bestVar = 0;

	double cumulMean = 0;
	double totalMean = 0;

	double cumulProb = 0;
	Vector<double> prob(255);

	//probabilities for each colour
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			++prob[(int)mtx(x,y)];
	prob /= (double)(mtx.height() * mtx.width());

	//total mean calculation
	for(uint i=0; i<prob.size(); ++i)
		totalMean += (i+1) * prob[i];

	//find best threshold
	for(double t=0; t<255; ++t)
	{
		cumulProb += prob[(uint)t];
		cumulMean += (t+1) * prob[(int)t];

		curVar = totalMean * cumulProb - cumulMean;
		curVar *= totalMean * cumulProb - cumulMean;
		curVar /= cumulProb * (1 - cumulProb);
		curVar = (cumulProb==0)? 0 : curVar;		//fix division by zero

		if(curVar > bestVar)
		{
			bestVar = curVar;
			threshold = t;
		}
	}

	NuLLProcessing::threshold(mtx, dst, threshold, gmin, gmax);
}

//simple gamma correction
template <typename T> void NuLLProcessing::gammaCorrection(const Matrix<T>& mtx, Matrix<T>& dst, double gamma)
{
	T maxVal = NuLLTools::maxValue(mtx);
	T gammaInv = 1.0 / gamma;

	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = maxVal * pow((mtx(x,y) / maxVal), gammaInv);
}

//contrast enhancement
template <typename T> void NuLLProcessing::logDynamicCompression(const Matrix<T>& mtx, Matrix<T>& dst, double c)
{
	if(!c)
	{
		c = 255;
		c /= log(1 + NuLLTools::maxValue(mtx));
	}

	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = c * log(1 + mtx(x,y));
}

//rescales all values to fit given range
template <typename T> void NuLLProcessing::affineRescale(const Matrix<T>& mtx, Matrix<T>& dst, double minVal, double maxVal)
{
	T val;
	T imgMax = NuLLTools::maxValue(mtx);
	T imgMin = NuLLTools::minValue(mtx);
	T imgDiff = (imgMax == imgMin)? imgMin : (imgMax - imgMin);		//if imgMax==imgMin, we get div by zero
	T valDiff = maxVal - minVal;

	for(uint y=0; y<mtx.height(); ++y)
	{
		for(uint x=0; x<mtx.width(); ++x)
		{
			val = (mtx(x,y) - imgMin) * (valDiff);
			val /= imgDiff;
			val += minVal;
			dst(x,y) = val;
		}
	}
}

//simple affine transformation
template <typename T> void NuLLProcessing::affineTransform(const Matrix<T>& mtx, Matrix<T>& dst, double a, double b)
{
	for(uint y=0; y<mtx.height(); ++y)
		for(uint x=0; x<mtx.width(); ++x)
			dst(x,y) = a * mtx(x,y) + b;
}

//create noise mask with uniform noise
template <typename T> void NuLLProcessing::uniformNoiseMask(Matrix<T>& dst, T a, T b, double percentage)
{

}

//create noise mask with gaussian noise
template <typename T> void NuLLProcessing::gaussianNoiseMask(Matrix<T>& dst, double sigma, double my, double percentage)
{

}

//create noise mask with salt and pepper noise
template <typename T> void NuLLProcessing::impulseNoiseMask(Matrix<T>& dst, double percentage)
{

}

#endif
