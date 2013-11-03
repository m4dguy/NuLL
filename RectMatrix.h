#ifndef RECTMATRIX_H
#define RECTMATRIX_H

#include "Matrix.h"

/*
 * quadratic matrix
 * plain and simple; works fine
 */

template <typename T> class RectMatrix : public Matrix<T>
{
    public:
        RectMatrix(size_t dimension)
        {
            this->_width = this->_height = dimension;
            this->_memory = dimension * dimension * sizeof(T);
            this->_entries = (T*) calloc(dimension * dimension, sizeof(T));
        };

        RectMatrix(size_t width, size_t height)
        {
            this->_width = width;
            this->_height = height;
            this->_memory = width * height * sizeof(T);
            this->_entries = (T*) calloc(width * height, sizeof(T));
        };

        RectMatrix(const RectMatrix& other)
        {
            this->_width = other._width;
            this->_height = other._height;
			this->_memory = other._memory;
            this->_entries = (T*) malloc(this->_memory);

            memcpy(this->_entries, other._entries, other._memory);
        };

        virtual ~RectMatrix(){free(this->_entries);};

        virtual T& get(uint x, uint y)
        {
            return this->_entries[y*this->_width + x];
        };

        virtual const T& get(uint x, uint y) const
        {
            return this->_entries[y*this->_width + x];
        };

        virtual void resize(size_t width, size_t height)
        {
            this->_width = width;
            this->_height = height;
            this->_memory = width * height * sizeof(T);

            this->_entries = (T*) realloc(this->_entries, this->_memory);
        };

        RectMatrix<T> operator*(const Matrix<T>& mtx)
        {
            T entry;
            size_t width = mtx.width();
            size_t height = mtx.height();
            RectMatrix<T> res(width, height);

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
    private:
};

typedef RectMatrix<double> RMatrixD;
typedef RectMatrix<float> RMatrixF;

#endif // RECTMATRIX_H
