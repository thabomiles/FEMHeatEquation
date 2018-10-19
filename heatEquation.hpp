#ifndef HEATEQUATIONHEADERREF
#define HEATEQUATIONHEADERREF
#include <iostream>
#include <vector>
#include <array>

class HeatEquation
{
public:
    void SetSpaceTimeMesh( std::vector<double> currentSpaceNodes, std::vector<double> currentTimeNodes );

    void SetVariables( int n, int m, double T, double a );

    void SetSystem( std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag,
        std::vector<double> f, std::vector<double> x );


private:
    void TridiagonalMatrixSolver( int n,
         std::vector<double> Diagonal, std::vector<double> LowerDiag,
         std::vector<double> UpperDiag, std::vector<double> f,
         std::vector<double> &x );

//    void TimeStepper();


    int mn, mm;

    double mT, ma;

    std::vector<double> mCurrentDiagonal;
    std::vector<double> mCurrentLowerDiag;
    std::vector<double> mCurrentUpperDiag;
    std::vector<double> mCurrentf;
    std::vector<double> mCurrentSolution;

    std::vector<double> mCurrentSpaceNodes;
    std::vector<double> mCurrentTimeNodes;




};


#endif // HEATEQUATIONHEADERREF
