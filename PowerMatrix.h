#ifndef POWERMATRIX_H
#define POWERMATRIX_H

#include "Matrix.h"

/*
 * matrix of dimension 2^n
 * fast access to elements due to bit operations
 * optimal memory usage for matrizes of size 2^n
 * for other dimensions, memory waste due to offcut
 * commonly used in different NuLL modules
 * use it like any normal matrix
 */

template <typename T> class PowerMatrix : public Matrix<T>
{
    public:
        PowerMatrix(size_t dimension)
        {
            this->_width = this->_height = dimension;

            this->_power = 1;
            uint trueSize = 1;
            while((trueSize<<=1) < dimension)
                ++(this->_power);

            this->_memory = trueSize * this->_height * sizeof(T);
            this->_entries = (T*) calloc(trueSize * this->_height, sizeof(T));
        };


        PowerMatrix(size_t width, size_t height)
        {
            this->_width = width;
            this->_height = height;

            this->_power = 1;
			uint trueSize = 1;
			while((trueSize<<=1) < width)
                ++(this->_power);

            this->_memory = trueSize * this->_height * sizeof(T);
            this->_entries = (T*) calloc(trueSize * this->_height, sizeof(T));
        };

        PowerMatrix(const PowerMatrix& other)
        {
            this->_width = other._width;
            this->_height = other._height;
			this->_memory = other._memory;
			this->_power = other._power;

            this->_entries = (T*) malloc(this->_memory);
            memcpy(this->_entries, other._entries, other._memory);
        };

        virtual ~PowerMatrix(){free(this->_entries);};

        virtual T& get(uint x, uint y)
        {
            return this->_entries[(y<<this->_power)|x];
        };

        virtual const T& get(uint x, uint y) const
        {
            return this->_entries[(y<<this->_power)|x];
        };

        virtual void resize(size_t width, size_t height)
        {
            this->_width = width;
            this->_height = height;
            this->_memory = width * height * sizeof(T);

            this->_power = 1;
            uint counter=1;
            while((counter<<1) < width)
                ++(this->_power);

            this->_entries = (T*) realloc(this->_entries, this->_memory);
        };

        PowerMatrix<T> operator*(const Matrix<T>& mtx)
        {
            T entry;
            size_t width = mtx.width();
            size_t height = mtx.height();
            PowerMatrix<T> res(width, height);

            for(uint j=0; j<height; ++j)
            {
                for(uint i=0; i<width; ++i)
                {
                    entry = 0;
                    for(uint k=0; k<height; ++k)
                    {
                        entry += this->get(k,i) * mtx(j,k);
                    }
                    res(j,i) = entry;
                }
            }
            return res;
        };

    protected:
        uint _power;

    private:
};

typedef PowerMatrix<double> PMatrixD;
typedef PowerMatrix<float> PMatrixF;

#endif // POWERMATRIX_H
