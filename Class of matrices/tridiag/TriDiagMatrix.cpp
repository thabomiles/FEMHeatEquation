#include <iostream>
#include <vector>
#include <array>
#include "TriDiagMatrix.hpp"

void TriDiagMatrix::SetMatrix( int n, std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag)
    {
        mpn = n;
        mpDiagonal = Diagonal;
        mpLowerDiag = LowerDiag;
        mpUpperDiag = UpperDiag;
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


//    void TridiagonalMatrixSolver( std::vector<double> f,
//         std::vector<double> &x );

