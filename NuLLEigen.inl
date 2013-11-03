#ifndef NULLEIGEN_INL
#define NULLEIGEN_INL

#include "NuLLEigen.h"

//calculates largest eigenvector and corresponding eigenvalue of a matrix
template<typename T> inline void NuLLEigen::powerIteration(const Matrix<T>& mtx, Vector<T>& dstVec, T& dstVal, T threshold)
{
    size_t dim = mtx.height();
    T lenOld, diff;
    dstVec[0] = 1;  //init; unnecessary?

    do
    {
        dstVec = mtx * dstVec;
        lenOld = dstVal;
        dstVal = NuLLMath::euclideanNorm(dstVec);

        dstVec /= dstVal;
        diff = (dstVal > lenOld)? (dstVal - lenOld) : (lenOld - dstVal);
    }
    while(diff > threshold);
}

//calculates all eigenvalues of a matrix
//uses QR decomposition and stores final eigenvalues in vector
//computationally expensive; do not use, if only the largest eigenvalue is needed!
template<typename T> inline void NuLLEigen::QRAlgorithm(const Matrix<T>& A, RectMatrix<T>& dstVecs, Vector<T>& dstVals, T threshold = 0.00001)
{
    size_t dim = A.height();

    T norm, normOld, diff;
    norm = normOld = diff = 0;

    RectMatrix<T> Ai(dim, dim);
    RectMatrix<T> Q(dim, dim);
    RectMatrix<T> R(dim, dim);
    NuLLTools::copyMatrix(A, Ai);

    do
    {
        normOld = norm;
        NuLLDecomposition::QRDecomposition(Ai, Q, R);
        Ai = R * Q;                                         //alternatively use:        NuLLMath::transpose(Q) * Ai * Q;
        //norm = NuLLMath::frobeniusNorm(Ai);

        norm = 0;
        for(uint i=0; i<dim; ++i)
            for(uint j=i; j<dim; ++j)
                norm += Ai(i,j);

        diff = (norm > normOld)? (norm - normOld) : (normOld - norm);
    }
    while(diff > threshold);

    for(uint i=0; i<dim; ++i)
        dstVals[i] = Ai(i,i);

    dstVecs = Q;
}

#endif
