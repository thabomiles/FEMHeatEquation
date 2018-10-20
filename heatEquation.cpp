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
        mCurrentTimeNodes = currentTimeNodes;

    }

    void HeatEquation::MatrixVectorMultiplier ( int n, std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag,
        std::vector<double> f, std::vector<double> &ProductVector )
    {
        ProductVector.assign( n, 0 );

        ProductVector.at(0) = Diagonal.at(0)*f.at(0)+UpperDiag.at(0)*f.at(1);
        ProductVector.at(n-1) = LowerDiag.at(n-1)*f.at(n-2)+Diagonal.at(n-1)*f.at(n-1);

            for (int j = 1; j<n-1; j++)
            {
                ProductVector.at(j) = LowerDiag.at(j)*f.at(j)+Diagonal.at(j)*f.at(j)+UpperDiag.at(j)*f.at(j+1);
            }

    }


//void HeatEquation::TimeStepper
