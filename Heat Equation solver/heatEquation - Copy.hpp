#ifndef HEATEQUATIONHEADERREF
#define HEATEQUATIONHEADERREF
#include <iostream>
#include <vector>
#include <array>

class HeatEquation
{
public:

    void Solve();

    void SetSpaceTimeMesh( std::vector<double> currentSpaceNodes, std::vector<double> currentTimeNodes );

    void SetVariables( int n, int m, double T, double a );

    void SetSystem( std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag,
        std::vector<double> f, std::vector<double> x );

    void MatrixVectorMultiplier ( int n, std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag,
        std::vector<double> f, std::vector<double> &ProductVector);

    void TridiagonalMatrixSolver( int n,
         std::vector<double> Diagonal, std::vector<double> LowerDiag,
         std::vector<double> UpperDiag, std::vector<double> f,
         std::vector<double> &x );

    void BuildMassMatrix ( int n, std::vector<double> SpaceMesh,
         std::vector<double> &MassDiagonal, std::vector<double> &MassLowerDiag,
         std::vector<double> &MassUpperDiag  );

    void BuildStiffnessMatrix ( int n, double a, std::vector<double> SpaceMesh,
         std::vector<double> &StiffnessDiagonal, std::vector<double> &StiffnessLowerDiag,
         std::vector<double> &StiffnessUpperDiag  );

private:

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
