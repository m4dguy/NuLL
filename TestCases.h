#include "SymMatrix.h"
#include "PowerMatrix.h"
#include "RectMatrix.h"

#include "NuLLMath.inl"
#include "NuLLEigen.inl"
#include "NuLLDecomposition.inl"
#include "NuLLSolve.inl"
#include "NuLLProcessing.inl"
#include "NuLLTools.inl"
#include "Vector.h"
#include "Utils.h"


//#define NDEBUG
#include <assert.h>

namespace TestCases
{
    //works sofar
    void testMathBaseRect()
    {
        RMatrixD I(3,3);
        I(0,0) = I(1,1) = I(2,2) = 1;

        RMatrixD A1(3,3);
        A1(0,0)=1;       A1(1,0)=2;       A1(2,0)=3;
        A1(0,1)=4;       A1(1,1)=5;       A1(2,1)=6;
        A1(0,2)=7;       A1(1,2)=8;       A1(2,2)=9;

        RMatrixD A2(3,3);
        A2(0,0)=3;       A2(1,0)=3;       A2(2,0)=3;
        A2(0,1)=8;       A2(1,1)=2;       A2(2,1)=8;
        A2(0,2)=0;       A2(1,2)=5;       A2(2,2)=4;

        VectorD V1(3);
        V1[0]=-1;        V1[1]=3;      V1[2]=0;        V1[3]=7;

        (A1 * I).print();
        (I * A1).print();
        (A1 * A1).print();
        (A1 * A2).print();
        (A2 * A1).print();
        //(A2 * V1).print();
    }

    //workz
    void testMathBaseRect2()
    {
        RMatrixD A1(3,3);
        A1(0,0)=1;       A1(1,0)=2;       A1(2,0)=3;
        A1(0,1)=4;       A1(1,1)=5;       A1(2,1)=6;
        A1(0,2)=7;       A1(1,2)=8;       A1(2,2)=9;
        A1.print();

        A1 *= 2.0;
        A1.print();
    }

    //fine
    void testMathBaseSym()
    {
        SMatrixD M1(3);
        M1(0,0)=1;
        M1(0,1)=4;      M1(1,1)=2;
        M1(0,2)=6;      M1(1,2)=5;        M1(2,2)=3;
        M1.print();

        SMatrixD I(3);
        I(0,0)=I(1,1)=I(2,2)=1.0;
        I.print();

        (I * M1).print();
        (M1 * I).print();
        (M1 * M1).print();
    }

    //works
    void testMathBaseSym2()
    {
        SMatrixD M1(3);
        M1(0,0)=1;
        M1(0,1)=4;      M1(1,1)=2;
        M1(0,2)=6;      M1(1,2)=5;        M1(2,2)=3;
        M1.print();

        RMatrixD M2(3);
        NuLLTools::copyMatrix(M1, M2);

        M2 *= 2.0;
        M2.print();

        M1 *= 2.0;
        M1.print();
    }

    //works fine
    void testTools()
    {
        RMatrixD A1(3,3);
        A1(0,0)=1;       A1(1,0)=2;       A1(2,0)=3;
        A1(0,1)=4;       A1(1,1)=5;       A1(2,1)=6;
        A1(0,2)=7;       A1(1,2)=8;       A1(2,2)=9;
        A1.print();

        SMatrixD M1(3);
        M1(0,0)=1;
        M1(0,1)=4;      M1(1,1)=2;
        M1(0,2)=6;      M1(1,2)=5;        M1(2,2)=3;
        M1.print();

        //VectorD V = NuLLTools::getColumn(A1, 2);
        //V.print();

        NuLLTools::copyMatrix(M1, A1);
        A1.print();
    }

    //test success
    void testQRDecomposition()
    {
        RMatrixD A6(4,4);
        A6(0,0)=-1;      A6(1,0)=1;      A6(2,0)=-2;      A6(3,0)=3;
        A6(0,1)=-1;      A6(1,1)=4;      A6(2,1)=-1;      A6(3,1)=6;
        A6(0,2)=1;       A6(1,2)=-4;     A6(2,2)=-1;      A6(3,2)=-4;
        A6(0,3)=2;       A6(1,3)=0;      A6(2,3)=5;       A6(3,3)=-4;
        std::cout << "original:" << std::endl << A6 << std::endl << std::endl;

        RMatrixD Q(A6.height());
        RMatrixD R(A6.height());
        NuLLDecomposition::QRDecomposition(A6, Q, R);

        std::cout << "product of Q and R:" << std::endl << (Q * R) << std::endl << std::endl;
    }

    //not fully working; trouble with cases where (0,0) != 1
    void testLUDecomposition()
    {
        PMatrixD mtx3(4,4);
        mtx3(0,0)=1;    mtx3(1,0)=-1;    mtx3(2,0)=-1;   mtx3(3,0)=1;
        mtx3(0,1)=-1;   mtx3(1,1)=-2;    mtx3(2,1)=0;    mtx3(3,1)=-2;
        mtx3(0,2)=-2;   mtx3(1,2)=-4;    mtx3(2,2)=2;    mtx3(3,2)=-2;
        mtx3(0,3)=1;    mtx3(1,3)=-7;    mtx3(2,3)=1;    mtx3(3,3)=1;

        /*PMatrixD mtx3(3,3);
        mtx3(0,0)=4;    mtx3(1,0)=4;     mtx3(2,0)=3;
        mtx3(0,1)=2;    mtx3(1,1)=4;    mtx3(2,1)=3;
        mtx3(0,2)=4;   mtx3(1,2)=3;     mtx3(2,2)=1;*/


        std::cout << "original:" << std::endl << mtx3 << std::endl << std::endl;

        PMatrixD L(mtx3.width(), mtx3.height());
        PMatrixD U(mtx3.width(), mtx3.height());
        NuLLDecomposition::LUDecomposition(mtx3, L, U);

        L.print();
        U.print();

        std::cout << "product of L and U:" << std::endl << (L * U) << std::endl << std::endl;
    }

    //fully working!
    void testCholeskyDecomposition()
    {
        SMatrixD A(4);
        A(0,0)=1;
        A(0,1)=5;      A(1,1)=29;
        A(0,2)=2;      A(1,2)=16;     A(2,2)=22;
        A(0,3)=0;      A(1,3)=2;      A(2,3)=15;     A(3,3)=33;
        std::cout << "original:" << std::endl << A << std::endl << std::endl;

        SMatrixD chol(A.height());
        NuLLDecomposition::choleskyDecomposition(A, chol);
        std::cout << "decomposed:" << std::endl << chol << std::endl << std::endl;
    }

    void testEigen()
    {
        RMatrixD mtx2(3,3);
        /*mtx2(0,0)=1;    mtx2(1,0)=1;    mtx2(2,0)=2;
        mtx2(0,1)=1;    mtx2(1,1)=2;    mtx2(2,1)=1;
        mtx2(0,2)=2;    mtx2(1,2)=1;    mtx2(2,2)=1;*/
        //eigenvalues = {4,-1,1}
        //eigenvectors = {{1,1,1},{-1,0,1},{1,-2,1}}

        mtx2(0,0)=2;    mtx2(1,0)=0;    mtx2(2,0)=0;
        mtx2(0,1)=0;    mtx2(1,1)=3;    mtx2(2,1)=4;
        mtx2(0,2)=0;    mtx2(1,2)=4;    mtx2(2,2)=9;

        VectorD vals(mtx2.height());
        RMatrixD vecs(mtx2.height());
        NuLLEigen::QRAlgorithm(mtx2, vecs, vals);

        std::cout << "vals:" << std::endl << vals << std::endl << std::endl;
        std::cout << "vecs:" << std::endl << vecs << std::endl << std::endl;
    }

    void testSolve()
    {
        /*

        SMatrixD mtx(6);
        mtx(0,0)=1;
        mtx(1,0)=2;     mtx(1,1)=5;
        mtx(2,0)=0;     mtx(2,1)=-1;    mtx(2,2)=2;
        mtx(3,0)=1;     mtx(3,1)=1;     mtx(3,2)=1;     mtx(3,3)=4;
        mtx(4,0)=1;     mtx(4,1)=2;     mtx(4,2)=-2;    mtx(4,3)=1;     mtx(4,4)=6;
        mtx(5,0)=-1;    mtx(5,1)=-2;    mtx(5,2)=1;     mtx(5,3)=1;     mtx(5,4)=-1;    mtx(5,5)=9;
        //mtx.print();

        RMatrixD mtx2(3);
        mtx2(0,0)=1;    mtx2(1,0)=2;    mtx2(2,0)=3;
        mtx2(0,1)=1;    mtx2(1,1)=1;    mtx2(2,1)=1;
        mtx2(0,2)=3;    mtx2(1,2)=3;    mtx2(2,2)=1;

        RMatrixD A2(4);
        A2(0,0)=1;      A2(0,1)=-1;     A2(0,2)=-1;     A2(0,3)=1;
        A2(1,0)=-1;     A2(1,1)=-2;     A2(1,2)=0;      A2(1,3)=-2;
        A2(2,0)=-2;     A2(2,1)=-4;     A2(2,2)=2;      A2(2,3)=-2;
        A2(3,0)=1;      A2(3,1)=-7;     A2(3,2)=1;      A2(3,3)=1;

        RMatrixD A3(5);
        A3(0,0)=5;      A3(0,1)=1;      A3(0,2)=0;      A3(0,3)=0;      A3(0,4)=0;
        A3(1,0)=1;      A3(1,1)=4;      A3(1,2)=1;      A3(1,3)=0;      A3(1,4)=0;
        A3(2,0)=0;      A3(2,1)=1;      A3(2,2)=6;      A3(2,3)=1;      A3(2,4)=0;
        A3(3,0)=0;      A3(3,1)=0;      A3(3,2)=1;      A3(3,3)=3;      A3(3,4)=1;
        A3(4,0)=0;      A3(4,1)=0;      A3(4,2)=0;      A3(4,3)=1;      A3(4,4)=3;

        RMatrixD A4(5);
        A4(0,0)=1;      A4(1,0)=0;      A4(2,0)=0;      A4(3,0)=0;      A4(4,0)=0;
        A4(0,1)=0;      A4(1,1)=1;      A4(2,1)=0;      A4(3,1)=0;      A4(4,1)=0;
        A4(0,2)=0;      A4(1,2)=0;      A4(2,2)=1;      A4(3,2)=0;      A4(4,2)=0;
        A4(0,3)=0;      A4(1,3)=0;      A4(2,3)=0;      A4(3,3)=1;      A4(4,3)=0;
        A4(0,4)=0;      A4(1,4)=0;      A4(2,4)=0;      A4(3,4)=0;      A4(4,4)=1;

        RMatrixD A5(4);
        A5(0,0)=20;      A5(1,0)=-113;   A5(2,0)=62;      A5(3,0)=9;
        A5(0,1)=-6;      A5(1,1)=32;     A5(2,1)=-17;     A5(3,1)=-4;
        A5(0,2)=-34;     A5(1,2)=182;    A5(2,2)=-97;     A5(3,2)=-18;
        A5(0,3)=92;      A5(1,3)=-481;   A5(2,3)=253;     A5(3,3)=43;

        VectorD vec(5);
        vec[0]=1;
        vec[1]=2;
        vec[2]=-1;
        vec[3]=9;
        vec[4]=9;
        //vec[5]=3;

        VectorD vec2(4);
        vec2[0]=-4;
        vec2[1]=-30;
        vec2[2]=-5;
        vec2[3]=67;


        //VectorD sol = NuLLSolve::conjugateGradient(mtx, vec);
        //VectorD sol1 = NuLLSolve::conjugateGradient(A4, vec);
        //VectorD sol2 = NuLLSolve::ThomasAlgorithm(A4, vec);
        //VectorD sol3 = NuLLSolve::solveLU(A5, vec2);

        */
    }


}
