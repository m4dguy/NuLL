#include "FLTKInterface.h"

//unchecked; probably troublesome/ error in implementation
/*void FLTKInterface::imageToPixelbuffer(const Fl_Image& src, uchar* dst)
{
	const size_t count = src.w() * src.h() * src.d();
	const char* data = src.data()[0];
	//std::memcpy(dst, data, sizeof(uchar));
	for(uint i=0; i<count; ++i)
		dst[i] = (uchar)data[i]; 
}*/

//assume buffer depth = 3
//fix for case of depth = 1
//assuming pixel buffer (called dst) exists
template<typename T> void FLTKInterface::matrixToPixelbuffer(const Matrix<T>& src, uchar* dst)
{
	//assume depth = 3
	uint index;				//buffer index
	const uint depth = 3;
	size_t height = src.height();
	size_t width = src.width();

	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			index = (y * width * depth) + (x * depth);
			dst[index] = (uchar) src(x,y);
			dst[index+1] = (uchar) src(x,y);
			dst[index+2] = (uchar) src(x,y);
		}
	}
}

//fix for case of depth = 1
template<typename T> void FLTKInterface::pixelbufferToMatrix(const uchar* src, Matrix<T>& dst)
{
	T val;
	uint index;				//buffer index
	const uint depth = 3;
	size_t height = dst.height();
	size_t width = dst.width();
	

	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			index = (y * width + x) * depth;
			val = (T)((src[index] + src[index+1] + src[index+2]) / 3.0);
			dst(x,y) = val;
		}
	}
}

//fix for case of depth = 1
template<typename T> void FLTKInterface::pixelbufferToMatrix(const uchar* src, Matrix<T>& dstR, Matrix<T>& dstG, Matrix<T>& dstB)
{
	uint index;				//buffer index
	size_t height = dstR.height();
	size_t width = dstR.width();

	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			index = (y * width + x) * depth;
			dstR(x,y) = (T)src[index];
			dstG(x,y) = (T)src[index+1];
			dstB(x,y) = (T)src[index+2];
		}
	}
}

template<typename T> void FLTKInterface::imageToMatrix(const Fl_Image& src, Matrix<T>& dst)
{
	uint index;				//buffer index
	const uint depth = src.d();
	uint height = src.h();
	uint width = src.w();
	
	T val;
	const char* buf = src.data()[0];
	uchar r,g,b;

	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			index = (y * width + x) * depth;
			switch(src.count())
			{
				case 1:
				{
					switch(depth)
					{
						case 1:
						{
							val = (T)buf[index];
							break;
						}
						case 3:
						{
							r = (uchar) buf[index];
							g = (uchar) buf[index+1];
							b = (uchar) buf[index+2];
							val = (T)((r + g + b) / 3.0);
							break;
						}
						default:
							printf("image depth not supported: chars=%d\n", depth);
							exit(1);
					}
					break;
				}
				default:
					printf("Not supported: count=%d\n", src.count());
					exit(1);
			}
 			dst(x,y) = val;
		}
	}
}

template<typename T> void FLTKInterface::imageToMatrix(const Fl_Image& src, Matrix<T>& dstR, Matrix<T>& dstG, Matrix<T>& dstB)
{
	uint index;				//buffer index
	const uint depth = src.d();
	uint height = src.height();
	uint width = src.width();

	const char* buf = src.data()[0];
	uchar r,g,b;

	for(uint y=0; y<height; ++y)
	{
		for(uint x=0; x<width; ++x)
		{
			index = (y * width * depth) + (x * depth);
			switch(src.count())
			{
				case 1:
				{
					switch(src.d())
					{
						case 1:
						{
							r = g = b = (uchar) buf[index];
							break;
						}
						case 3:
						{
							r = (uchar) buf[index];
							g = (uchar) buf[index+1];
							b = (uchar) buf[index+2];
							break;
						}
						default:
							printf("Not supported: chars=%d\n", src.d());
							exit(1);
					}
					break;
				}
				default:
					printf("Not supported: count=%d\n", src.count());
					exit(1);
			}
			
			dstR(x,y) = (T)r;
			dstG(x,y) = (T)g;
			dstB(x,y) = (T)b;
		}
	}
}
