#ifndef NULLDECOMPOSITION_INL
#define NULLDECOMPOSITION_INL

#include "NuLLDecomposition.h"

//works only for symmetric positive definite matrizes
//useful for linear least squares, monte carlo, ...
//for solving linear systems, use conjugate gradients instead!
//remember: cholesky matrices are triangular!
template <typename T> void NuLLDecomposition::choleskyDecomposition(const Matrix<T>& mtx, Matrix<T>& dst)
{
    size_t dim = mtx.height();
    T sum, arg;

    //first column
    dst(0,0) = sqrt(mtx(0,0));
    for(uint i=1; i<dim; ++i)
    {
        arg = mtx(i,0);
        arg /= dst(0,0);
        dst(i,0) = arg;
    }

    //calclation row by row
    for(uint i=1; i<dim; ++i)
    {
        //diagonal entry
        arg = mtx(i,i);
        sum = 0;
        for(uint j=0; j<=i-1; ++j)
        {
            sum += dst(i,j) * dst(i,j);
        }
        arg -= sum;
        dst(i,i) = sqrt(arg);

        //off-diagonals
        for(uint j=i+1; j<dim; ++j)
        {
            arg = mtx(i,j);

            sum = 0;
            for(uint k=0; k<=j-1; ++k)
            {
                sum += dst(j,k) * dst(i,k);
            }
            arg -= sum;
            arg /= dst(i,i);
            dst(j,i) = arg;
        }
    }
}

//TODO: debug! dont use temporary matrix!
//trouble with cases where (0,0) != 1
//decomposition applicable to all invertible matrizes
//L and U are triangular matrizes (L is lower left triangular, U is upper right triangular)
//use for solving linear equations, if CG is not applicable!
//used for matrix inversion and determinant
template <typename T> void NuLLDecomposition::LUDecomposition(const Matrix<T>& mtx, Matrix<T>& dstL, Matrix<T>& dstU)
{
    size_t dim = mtx.width();
    PowerMatrix<T> tmp(dim, dim);
    NuLLTools::copyMatrix(mtx, tmp);

    for(uint y=1; y<dim; ++y)
        tmp(0,y) /= tmp(0,0);

    for(uint y=0; y<dim-1; ++y)
    {
        for(uint x=y+1; x<dim; ++x)
        {
            tmp(y,x) /= tmp(y,y);
            for(uint k=y+1; k<dim; ++k)
            {
                tmp(k,x) -= (tmp(k,y) * tmp(y,x));
            }
        }
    }

    dstL.fill();
    dstU.fill();
    for(uint y=0; y<dim; ++y)
    {
        dstL(y,y) = 1.0;
        for(uint x=0; x<dim; ++x)
        {
            if(x<y)
                dstL(x,y) = tmp(x,y);
            else
                dstU(x,y) = tmp(x,y);
        }
    }
}

//optimize: avoid binary operators!
//internally, PowerMatrices are used to speed up calculations
//return orthogonal basis Q and upper triangular matrix R
//useful for QR algorithm (eigenvalue calculation)
template <typename T> void NuLLDecomposition::QRDecomposition(const Matrix<T>& mtx, Matrix<T>& dstQ, Matrix<T>& dstR)
{
    T s = 0;
    T w = 0;
    T tmp = 0;
    size_t dim = mtx.height();

	Vector<T> h(dim);

	PowerMatrix<T> Q(dim);				//Q
	PowerMatrix<T> R(dim);				//R
    PowerMatrix<T> I(dim);				//unit matrix
    PowerMatrix<T> Hk(dim);				//householder matrix H of step k
    PowerMatrix<T> hh(dim);				//dyadic product of h
    NuLLTools::copyMatrix(mtx, R);

    //build unit matrix, initialize orthogonal matrix
    for(uint i=0; i<dim; ++i)
        I(i,i) = Q(i,i) = 1;

    for(uint k=0; k<dim-1; ++k)
    {
        //preparation for householder matrix
        //calculation of vector h
        s = 0;
        for(uint j=k; j<dim; ++j)
        {
            s += R(k,j) * R(k,j);
        }

        s = sqrt(s);
        if(R(k,k) >= 0)
        {
            s *= -1;
        }

        //building vector for householder matrix
        for(uint i=0; i<k; ++i)
        {
            h[i] = 0;
        }

        tmp = .5 * (1 - R(k,k) / s);
        h[k] = sqrt(tmp);

        tmp = 2 * h[k] * s;
        w = -1 / tmp;
        for(uint j=k+1; j<dim; ++j)
        {
            h[j] = w * R(k,j);
        }

        //calculation of householder matrix
        Hk = I;
        hh *= 2;
        Hk -= hh;
        R = Hk * R;
        Q = Q * Hk;
    }

	NuLLTools::copyMatrix(Q, dstQ);
	NuLLTools::copyMatrix(R, dstR);
}


#endif
