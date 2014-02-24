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
	//T col;

	const uint width = img.width();
	const uint height = img.height();

	for(uint i=0; i<circles.size(); ++i)
	{
		HoughCircle c = circles[i];
		x = c.x;
		y = c.y;
		rad = c.radius;
		radSqLow = (c.radius - .5) * (c.radius - .5);
		radSqUp = (c.radius + .5) * (c.radius + .5);
		//col = (T)c.score;

		for(int ry=-rad; ry<=rad; ++ry)
		{
			for(int rx=-rad; rx<=rad; ++rx)
			{
				radSq = (rx*rx + ry*ry);
				if((radSqLow <= radSq) && (radSq < radSqUp))
				{
					img(x+rx, y+ry) = 255;
				}
				
			}
		}
	}
}

void HoughTransform::mergeCircles(std::vector<HoughCircle>& circles, uint maxDist)
{
	uint dist;
	uint bSize;
	const uint maxDistSq = maxDist * maxDist;

	HoughCircle c;

	std::vector<uint> done;
	std::vector<HoughCircle> batch;
	std::vector<HoughCircle> res;
	done.resize(circles.size());
	batch.reserve(circles.size());
	res.reserve(circles.size());

	for(uint i=0; i<circles.size()-1; ++i)
	{
		for(uint j=i+1; j<circles.size(); ++j)
		{
			if(done[j])
				continue;

			dist = (circles[i].x - circles[j].x) * (circles[i].x - circles[j].x) + (circles[i].y - circles[j].y) * (circles[i].y - circles[j].y);
			if(dist <= maxDistSq)
			{
				done[i] = 1;
				batch.push_back(circles[j]);
			}
		}

		bSize = batch.size();
		if(bSize)
		{
			c = batch[0];
			for(uint b=0; b<bSize; ++b)
			{
				if(batch[b].score > c.score)
					c = batch[b];
			}
			res.push_back(c);
			batch.clear();
		}
	}
	res.swap(circles);
}

template<typename T> void HoughTransform::histogramInCircle(const Matrix<T>& mtx, const HoughCircle& circ, std::vector<uint> hist)
{
	uint x = circ.x;
	uint y = circ.y;
	int radius = circ.radius;
	int radSq = radius * radius;
	
	hist.clear();
	hist.reserve(256);
	hist.resize(256);

	for(int ry=-radius; ry<=radius; ++ry)
	{
		for(int rx=-radius; rx<=radius; ++rx)
		{
			if((rx*rx + ry*ry) <= radSq)
			{
				++hist[(int)mtx(x+rx, y+ry)];
			}
		}
	}

}
