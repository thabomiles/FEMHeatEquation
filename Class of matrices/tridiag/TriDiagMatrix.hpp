#ifndef TRIDIAGMATRIXHEADERREF
#define TRIDIAGMATRIXHEADERREF
#include <iostream>
#include <vector>
#include <array>

class TriDiagMatrix
{
public:

    void SetMatrix( int n, std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag);

    void MatrixVectorMultiplier ( std::vector<double> f, std::vector<double> &ProductVector);

    void PrintMatrix();

//    void TridiagonalMatrixSolver( std::vector<double> f,
//         std::vector<double> &x );

    void BuildMassMatrix ();

    void BuildStiffnessMatrix ();

private:

    int mpn;

    std::vector<double> mpDiagonal;
    std::vector<double> mpLowerDiag;
    std::vector<double> mpUpperDiag;

};


#endif // TRIDIAGMATRIXHEADERREF

