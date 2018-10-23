#include <iostream>
#include <cmath>
#include "heatEquation.hpp"
#include <vector>
#include <array>
#include <math.h>

void HeatEquation::Solve()
{

}

void HeatEquation::TridiagonalMatrixSolver( int n,
         std::vector<double> Diagonal, std::vector<double> LowerDiag,
         std::vector<double> UpperDiag, std::vector<double> f,
         std::vector<double> &x )
    {
                  //scaling matrix entries
    for(int i=1; i<n; i++)
    {
        Diagonal.at(i) = Diagonal.at(i)-(UpperDiag.at(i-1)*LowerDiag.at(i)/Diagonal.at(i-1));
        f.at(i) = f.at(i)-(f.at(i-1)*LowerDiag.at(i)/Diagonal.at(i-1));
    }

        //solving through elimination
     x.at(n-1)=f.at(n-1)/Diagonal.at(n-1);
     for(int i=n-2; i>=0; i=i-1)
     {
         x.at(i) = (f.at(i) - UpperDiag.at(i)*x.at(i+1))/(Diagonal.at(i));
     }
    }


void HeatEquation::SetSystem( std::vector<double> Diagonal,
            std::vector<double> LowerDiag, std::vector<double> UpperDiag,
            std::vector<double> f, std::vector<double> x )
    {
        mCurrentf = f;
        mCurrentSolution = x;

    }


//void HeatEquation::TimeStepper
