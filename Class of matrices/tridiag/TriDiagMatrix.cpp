#include <iostream>
#include <vector>
#include <array>
#include "TriDiagMatrix.hpp"

void TriDiagMatrix::SetMatrix(  std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag)
    {

        mpDiagonal = Diagonal;
        mpLowerDiag = LowerDiag;
        mpUpperDiag = UpperDiag;

        if ( (mpDiagonal.size()==mpLowerDiag.size())
            &&(mpDiagonal.size()==mpUpperDiag.size()))
        {
            mpn = mpDiagonal.size();
        }
        else
        {
            std::cout << ' your vectors are different sizes ';
        }

    }

void TriDiagMatrix::PrintMatrix()
{
        for (auto j: mpUpperDiag)
        std::cout << j << ' ';
        std::cout << " \n";

        for (auto j: mpDiagonal)
        std::cout << j << ' ';
        std::cout << " \n";

        for (auto j: mpLowerDiag)
        std::cout << j << ' ';
        std::cout << " \n";
}

void TriDiagMatrix::MatrixVectorMultiplier ( std::vector<double> f, std::vector<double> &ProductVector )
    {
        ProductVector.assign( mpn, 0 );

        ProductVector[0] = mpDiagonal[0]*f[0]+mpUpperDiag[0]*f[1];
        ProductVector[mpn-1] = mpLowerDiag[mpn-1]*f[mpn-2]+mpDiagonal[mpn-1]*f[mpn-1];

            for (int j = 1; j<mpn-1; j++)
            {
                ProductVector[j] = mpLowerDiag[j]*f[j-1]+mpDiagonal[j]*f[j]+mpUpperDiag[j+1]*f[j+1];
            }
    }


    void TriDiagMatrix::MatrixSolver( std::vector<double> f, std::vector<double> &x )
    {
        std::vector<double> cDiagonal = mpDiagonal;
        std::vector<double> cLowerDiag = mpLowerDiag;
        std::vector<double> cUpperDiag = mpUpperDiag;
         //scaling matrix entries
    for(int i=1; i<mpn; i++)
    {
        cDiagonal.at(i) = cDiagonal.at(i)-(cUpperDiag.at(i-1)*cLowerDiag.at(i)/cDiagonal.at(i-1));
        f.at(i) = f.at(i)-(f.at(i-1)*cLowerDiag.at(i)/cDiagonal.at(i-1));
    }

        //solving through elimination
     x.at(mpn-1)=f.at(mpn-1)/cDiagonal.at(mpn-1);
     for(int i=mpn-2; i>=0; i=i-1)
     {
         x.at(i) = (f.at(i) - cUpperDiag.at(i)*x.at(i+1))/(cDiagonal.at(i));
     }
    }

    void TriDiagMatrix::MultiplyByScalar( double k )
    {
        for (int i = 0; i<mpn; i++)
        {
            mpDiagonal.at(i) = k*mpDiagonal.at(i);
            mpLowerDiag.at(i) = k*mpLowerDiag.at(i);
            mpUpperDiag.at(i) = k*mpUpperDiag.at(i);
        }
    }

        //the sum initialises a new matrix
    void TriDiagMatrix::AddTwoMatrices ( TriDiagMatrix firstMatrix, TriDiagMatrix secondMatix )
    {
        if (firstMatrix.mpn==secondMatix.mpn)
        {
            mpn = secondMatix.mpn;
            mpDiagonal.clear();
            mpLowerDiag.clear();
            mpUpperDiag.clear();

            for (int i = 0; i<mpn; i++)
            {
                mpDiagonal.push_back( firstMatrix.mpDiagonal.at(i)+secondMatix.mpDiagonal.at(i) );
                mpLowerDiag.push_back(firstMatrix.mpUpperDiag.at(i)+secondMatix.mpUpperDiag.at(i) );
                mpUpperDiag.push_back( firstMatrix.mpLowerDiag.at(i)+secondMatix.mpLowerDiag.at(i) );
            }
        }
        else
        {
            std::cout << ' your vmatrices are different sizes ';
        }
    }

