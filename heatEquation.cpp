#include <iostream>
#include "heatEquation.hpp"
#include <vector>
#include <array>

void HeatEquation::TridiagonalMatrixSolver( int n,
         std::vector<double> Diagonal, std::vector<double> LowerDiag,
         std::vector<double> UpperDiag, std::vector<double> f,
         std::vector<double> &x )
    {
                  //scaling matrix entries
    for(int i=1; i<n; i++)
    {
        Diagonal[i] = Diagonal[i]-(UpperDiag[i-1]*LowerDiag[i]/Diagonal[i-1]);
        f[i] = f[i]-(f[i-1]*LowerDiag[i]/Diagonal[i-1]);
    }

        //solving through elimination
     x[n-1]=f[n-1]/Diagonal[n-1];
     for(int i=n-2; i>=0; i=i-1)
     {
         x[i] = (f[i] - UpperDiag[i]*x[i+1])/(Diagonal[i]);
     }
    }


    void HeatEquation::SetSystem( std::vector<double> Diagonal,
            std::vector<double> LowerDiag, std::vector<double> UpperDiag,
            std::vector<double> f, std::vector<double> x )
    {
        mCurrentDiagonal = Diagonal;
        mCurrentLowerDiag = LowerDiag;
        mCurrentUpperDiag = UpperDiag;
        mCurrentf = f;
        mCurrentSolution = x;

    }


    void HeatEquation::SetVariables( int n, int m, double T, double a )
    {
        mn = n;
        mm = m;
        mT = T;
        ma = a;

    }

    void HeatEquation::SetSpaceTimeMesh( std::vector<double> currentSpaceNodes,
                           std::vector<double> currentTimeNodes )
    {
        mCurrentSpaceNodes = currentSpaceNodes;
        mCurrentTimeNodes = currentTimeNodes:

    }



//void HeatEquation::TimeStepper
