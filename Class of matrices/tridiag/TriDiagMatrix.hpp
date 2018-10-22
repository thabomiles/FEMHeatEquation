#ifndef TRIDIAGMATRIXHEADERREF
#define TRIDIAGMATRIXHEADERREF
#include <iostream>
#include <vector>
#include <array>

class TriDiagMatrix
{
public:

    void SetMatrix( std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag);

    void MatrixVectorMultiplier ( std::vector<double> f, std::vector<double> &ProductVector);

    void PrintMatrix();

    void MatrixSolver( std::vector<double> f,
         std::vector<double> &x );

    void MultiplyByScalar (double k);

    void AddTwoMatrices ( TriDiagMatrix firstMatrix, TriDiagMatrix secondMatix );

protected:

    int mpn;

    std::vector<double> mpDiagonal;
    std::vector<double> mpLowerDiag;
    std::vector<double> mpUpperDiag;

};


#endif // TRIDIAGMATRIXHEADERREF

