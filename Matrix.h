#ifndef MATRIX_H
#define MATRIX_H

#include <string.h>

#include "Utils.h"
#include "Vector.h"

/*
 * matrix base class
 * needs overloaded get() for entries
 */

template <typename T> class Matrix
{
    public:
        Matrix(size_t dimension)
        {
            this->_width = this->_height = dimension;

            this->_power = 1;
            uint trueSize = 1;
            while((trueSize<<=1) < dimension)
                ++(this->_power);

            this->_memory = trueSize * this->_height * sizeof(T);
            this->_entries = (T*) calloc(trueSize * this->_height, sizeof(T));
        };


        Matrix(size_t width, size_t height)
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

        Matrix(const Matrix& other)
        {
            this->_width = other._width;
            this->_height = other._height;
			this->_memory = other._memory;
			this->_power = other._power;

            this->_entries = (T*) malloc(this->_memory);
            memcpy(this->_entries, other._entries, other._memory);
        };

        virtual ~Matrix()
        {
            free(this->_entries);
        };

        virtual T& get(uint x, uint y)
        {
            return this->_entries[(y<<this->_power)|x];
        };

        virtual const T& get(uint x, uint y) const
        {
            return this->_entries[(y<<this->_power)|x];
        };

        const T& getMirrored(int x, int y) const
        {
            int w = _width;
			int h = _height;
			x = abs(x);
            y = abs(y);
            x = (x>=w)? (w - (x - w + 1)) : x;
            y = (y>=h)? (h - (y - h + 1)) : y;
            return get(x,y);
        };

        T& getMirrored(int x, int y)
        {
            int w = _width;
			int h = _height;
			x = abs(x);
            y = abs(y);
            x = (x>=w)? (w - (x - w + 1)) : x;
            y = (y>=h)? (h - (y - h + 1)) : y;

            return get(x,y);
        };

        inline void set(uint x, uint y, T& val)
        {
            (*this)(x,y) = val;
        };

        inline size_t width() const
        {
            return _width;
        };

        inline size_t height() const
        {
            return _height;
        };

		void fill(T val = 0)
		{
			for(uint y=0; y<_height; ++y)
				for(uint x=0; x<_width; ++x)
					get(x,y) = val;
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

        inline T sum()
        {
            T res = 0;
            for(uint x=0; x<this->_width; ++x)
                for(uint y=0; y<this->_height; ++y)
                    res += this->get(x,y);

            return res;
        }

        inline T& operator()(uint x, uint y)
        {
            return get(x,y);
        }

        inline const T& operator()(uint x, uint y) const
        {
            return get(x,y);
        };

        const Matrix& operator=(const Matrix<T>& other)
        {
            if(this == &other)   // handle self assignment
                return *this;

            //reallocate if not enough memory
            if(other._memory > _memory)
                _entries = (T*) realloc(_entries, other._memory);

            _memory = other._memory;

            _width = other._width;
            _height = other._height;
            memcpy(_entries, other._entries, _memory);

            return *this;
        };

        Matrix<T>& operator+=(const Matrix<T>& A)
        {
            for(uint x=0; x<this->_width; ++x)
                for(uint y=0; y<this->_height; ++y)
                    this->get(x,y) += A(x,y);

            return *this;
        };

        Matrix<T>& operator-=(const Matrix<T>& A)
        {
            for(uint x=0; x<this->_width; ++x)
                for(uint y=0; y<this->_height; ++y)
                    this->get(x,y) -= A(x,y);

            return *this;
        };

        Matrix<T> operator*(const Matrix<T>& mtx)
        {
            T entry;
            size_t width = mtx.width();
            size_t height = mtx.height();
            Matrix<T> res(width, height);

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

        Vector<T> operator*(const Vector<T>& vec)
        {
            T entry;
            Vector<T> res(_width);

            for(uint x=0; x<_width; ++x)
            {
                entry = 0;
                for(uint y=0; y<_height; ++y)
                {
                    entry += vec[y] * get(x,y);
                }
                res[x] = entry;
            }

            return res;
        };

        Matrix<T>& operator*=(T scalar)
        {
            size_t dim = _memory / sizeof(T);
            for(uint i=0; i<dim; ++i)
                _entries[i] *= scalar;

            return *this;
        };

		//optimize!
        Matrix<T>& operator*=(const Matrix<T>& A)
        {
            T entry;
            size_t w = A.width();
            size_t h = A.height();
            Vector<T> row(w);

            for(uint y=0; y<h; ++y)
            {
                for(uint x=0; x<w; ++x)
                {
                    //save row; will be overwritten otherwise
                    for(uint k=0; k<w; ++k)
                        row[k] = this->get(k,y);

                    entry = 0;
                    for(uint k=0; k<h; ++k)
                        entry += row[k] * A(x,k);

                    this->get(y,x) = entry;
                }
            }
            return *this;
        };

        Matrix<T>& operator/=(T scalar)
        {
            size_t dim = _memory / sizeof(T);
            for(uint i=0; i<dim; ++i)
                _entries[i] /= scalar;

            return *this;
        };

        bool operator==(const Matrix<T>& other) const
        {
            if((other._width != this->_width) && (other._height != this->_height))
                return 0;

            bool eq = 1;
            for(uint x=0; x<_width; ++x)
                for(uint y=0; y<_height; ++y)
                    eq &= (other.get(x,y) == get(x,y));

            return eq;
        };

        bool operator!=(const Matrix<T>& other) const
        {
            return !(*this == other);
        };

		//output in csv format for file dumps
        friend std::ostream& operator<<(std::ostream& stream, const Matrix<T>& mtx)
        {
            for(uint y=0; y<mtx._height; ++y)
            {
                for(uint x=0; x<mtx._width; ++x)
                {
                    stream << mtx.get(x,y);

					if(mtx._width - x - 1)
						stream << ", ";
				}
                stream << std::endl;
            }
            return stream;
        };

        //use for debugging
        inline void print()
        {
            std::cout << (*this) << std::endl;
        };

    protected:
        uint _power;
        size_t _width;          //matrix dimensions
        size_t _height;

        size_t _memory;         //used for memory management and operators
        T* _entries;            //pointer to entries

    private:
};

typedef Matrix<double> MatrixD;
typedef Matrix<float> MatrixF;

#endif // MATRIX_H
