#ifndef SYMMATRIX_H
#define SYMMATRIX_H

#include <algorithm>

#include "RectMatrix.h"

/*
 * Symmetric matrix
 * slow access but good memory usage
 * only quadratic dimension possible
 *
 * fix: operator*= broken due to memory management
 *      progenitor of fixed version included. fix problem with use of template<template>
 *
 */

template <typename T> class SymMatrix : public Matrix<T>
{
    public:
        SymMatrix(size_t dim)
        {
            this->_width = this->_height = dim;

            this->_memory = (dim*(dim+1)>>1) * sizeof(T);
            this->_entries = (T*) calloc(dim*(dim+1)>>1, sizeof(T));
        };

        SymMatrix(const SymMatrix& other)
        {
            this->_width = this->_height = other._width;
            this->_memory = other._memory;
            this->_entries = (T*) malloc(other._memory);

            memcpy(this->_entries, other._entries, other._memory);
        };

        virtual ~SymMatrix(){free(this->_entries);};

        virtual T& get(uint x, uint y)
        {
            uint pos;
            if(x<y)
                std::swap(x,y);

            //umoptimized formula for index calculation:
            //pos = (x*(x+1)/2)+y
            pos = (x*(x+1)>>1) + y;
            return this->_entries[pos];
        };

        virtual const T& get(uint x, uint y) const
        {
            uint pos;
            if(x<y)
                std::swap(x,y);

            pos = (x*(x+1)>>1) + y;
            return this->_entries[pos];
        };

        //second parameter is ignored
        virtual void resize(size_t dim, size_t)
        {
            this->_width = this->_height = dim;
            this->_memory = (dim*(dim+1)>>1) * sizeof(T);

            this->_entries = (T*) realloc(this->_entries, this->_memory);
        };

        RectMatrix<T> operator*(const Matrix<T>& mtx)
        {
            T entry;
            size_t dim = mtx.width();
            RectMatrix<T> res(dim, dim);

            for(uint j=0; j<dim; ++j)
            {
                for(uint i=0; i<dim; ++i)
                {
                    entry = 0;
                    for(uint k=0; k<dim; ++k)
                    {
                        entry += this->get(k,i) * mtx(j,k);
                    }
                    res(j,i) = entry;
                }
            }
            return res;
        };

        //implement
        /*SymMatrix<T>& operator*=(const Matrix<T>& mtx)
        {
            T entry;
            size_t dim = mtx.width();
            Vector<T> row(dim);

            for(uint y=0; y<dim; ++y)
            {
                for(uint x=0; x<dim; ++x)
                {
                    //save row; will be overwritten otherwise
                    for(uint k=0; k<dim; ++k)
                        row[k] = this->get(k,y);

                    entry = 0;
                    for(uint k=0; k<dim; ++k)
                        entry += row[k] * mtx(x,k);

                    this->get(y,x) = entry;
                }
            }
            return *this;
        };*/

    protected:
    private:
};

typedef SymMatrix<double> SMatrixD;
typedef SymMatrix<float> SMatrixF;

#endif // SYMMATRIX_H
