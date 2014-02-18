#include "HoughTransform.h"

template<typename T> void HoughTransform::houghCircles(const Matrix<T>& img, std::vector<HoughCircle>& circles, double threshold, int minRad, int maxRad)
{
	HoughCircle circBest;
	uint scale;
	double score;
	double radSq, radSqLow, radSqUp;

	const uint width = img.width();
	const uint height = img.height();

	for(uint y=minRad; y<height-minRad; y+=2)
	{
		for(uint x=minRad; x<width-minRad; x+=2)
		{
			circBest.score = 0;
			//find circles of varying radius
			for(int rad=minRad; rad<=maxRad; ++rad)
			{
				if((x+rad >= width) || (y+rad >= height) || (x-rad <= 0) || (y-rad <= 0))
					break;

				radSqLow = (rad - .5) * (rad - .5);
				radSqUp = (rad + .5) * (rad + .5);
				score = 0;
				scale = 1;
				//scan area
				for(int ry=-rad; ry<=rad; ++ry)
				{
					for(int rx=-rad; rx<=rad; ++rx)
					{
						radSq = (rx*rx + ry*ry);
						if((radSqLow <= radSq) && (radSq < radSqUp))
						{
							++scale;
							score += img(x+rx, y+ry);
						}
					}
				}
				
				//check if circle found
				score /= scale;
				if(score > circBest.score)
				{
					circBest.x = x;
					circBest.y = y;
					circBest.radius = rad;
					circBest.score = score;
				}
			}

			if(circBest.score > threshold)
			{
				circles.push_back(circBest);
			}
		
		}
	}
}

template<typename T> void HoughTransform::plotCircles(const std::vector<HoughCircle>& circles, Matrix<T>& img)
{
	uint x, y;
	int rad, radSq;
	double radSqLow, radSqUp;
	T col;
	uint width = img.width();
	uint height = img.height();

	for(uint i=0; i<circles.size(); ++i)
	{
		HoughCircle c = circles[i];
		x = c.x;
		y = c.y;
		rad = c.radius;
		radSqLow = (c.radius - .5) * (c.radius - .5);
		radSqUp = (c.radius + .5) * (c.radius + .5);
		col = (T)c.score;

		for(int ry=-rad; ry<=rad; ++ry)
		{
			for(int rx=-rad; rx<=rad; ++rx)
			{
				radSq = (rx*rx + ry*ry);
				if((radSqLow <= radSq) && (radSq < radSqUp))
				{
					img(x+rx, y+ry) = col;
				}
				
			}
		}
	}
}
