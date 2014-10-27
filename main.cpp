#include "Matrix.h"
#include "Vector.h"
#include "Utils.h"


//#include "NuLLMath.inl"
//#include "NuLLEigen.inl"
//#include "NuLLDecomposition.inl"
//#include "NuLLSolve.inl"
//#include "NuLLProcessing.inl"
//#include "NuLLTools.inl"

//#include "TestCases.h"
#include "NuLLTransform.inl"


/*
 * NuLL: Numerical Lightweight Library
 * use for: solving linear systems, matrix decomposition,
 *          eigenvalue decomposition, data processing
 */

int main()
{
    MatrixD A(2,2);
    A(0,0)=2;   A(1,0)=1;
    A(0,1)=5;   A(1,1)=7;
    A.print();

    MatrixD B(2,2);
    MatrixD C(2,2);

    NuLLTransform::fourierTransform(A,B,C);
    B.print();

    //NuLLTransform::cosineTransform(B,A);
    //A.print();

    /*VectorD b(2);
    b[0]=11;
    b[1]=13;
    b.print();

    VectorD x(2);
    x[0]=1;
    x[1]=1;*/

    /*MatrixD A(3,3);
    A(0,0)=2;   A(1,0)=0;   A(2,0)=1;
    A(0,1)=0;   A(1,1)=1;   A(2,2)=0;
    A(0,2)=1;   A(1,2)=0;   A(2,2)=1;
    A.print();

    VectorD b(3);
    b[0]=0;
    b[1]=1;
    b[2]=1;
    b.print();

    VectorD x(3);
    x[0]=-1;
    x[1]=0;
    x[2]=0;*/

    //solution = 7.111, -3.222
    /*NuLLSolve::gaussSeidelSOR(A, b, x);
    x.print();

    MatrixD tmp(2,2);
    tmp = A;

    double lambda = 0;
    VectorD v(2);
    v[0] = v[1] = 4;
    NuLLEigen::powerIteration(A, v, lambda);
    NuLLTools::normalize(v);
    v.print();

    double lambdaInv = 1. / lambda;
    x = A * v;
    x *= lambdaInv;
    x.print();*/




    //QMatrixD m(1,1);
    /*PMatrixD m(4,4);
    m(0,0)=1;    m(1,0)=2;  m(2,0)=3;  m(3,0)=4;
    m(0,1)=5;    m(1,1)=6;  m(2,1)=7;  m(3,1)=8;
    m(0,2)=9;    m(1,2)=10; m(2,2)=11; m(3,2)=12;
    m(0,3)=13;   m(1,3)=14; m(2,3)=15; m(3,3)=16;
    m.print();*/

    /*MatrixD A(3,2);
    A(0,0)=1;   A(1,0)=2;   A(2,0)=5;
    A(0,1)=3;   A(1,1)=4;   A(2,1)=6;
    A.print();*/

    /*RMatrixD B(2,2);
    B(0,0)=7;   B(1,0)=8;
    B(1,0)=9;   B(1,1)=0;

    RMatrixD C(4,6);

    NuLLMath::kroneckerProduct(A,B,C);
    C.print();*/

    //TestCases::testMathBaseQuad();
    //TestCases::testMathBaseSym();

    //TestCases::testMathBaseQuad2();
    //TestCases::testMathBaseSym();

    //TestCases::testEigen();
    //TestCases::testTools();
    //TestCases::testLUDecomposition();

    /*Vector<double> vec1(4);
    Vector<double> vec2(4);*/

    /*MatrixD mtx(4,4);
    mtx(0,0)=1;  mtx(1,0)=0;  mtx(2,0)=2;  mtx(3,0)=-1;
    mtx(0,1)=3;  mtx(1,1)=0;  mtx(2,1)=0;  mtx(3,1)=5;
    mtx(0,2)=2;  mtx(1,2)=1;  mtx(2,2)=4;  mtx(3,2)=-3;
    mtx(0,3)=1;  mtx(1,3)=0;  mtx(2,3)=5;  mtx(3,3)=0;
    mtx.print();*/

    //vec2 = mtx * vec1;

    //NuLLSolve::conjugateGradient(mtx, vec1, vec2);

    //double det = NuLLMath::determinant(mtx);
    //std::cout << det << std::endl;
    //TestCases::testCholeskyDecomposition();
    //TestCases::testQRDecomposition();

    return 0;
}
