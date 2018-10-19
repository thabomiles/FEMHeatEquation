#include <iostream>
#include "heatEquation.hpp"

void HeatEquation::TridiagonalMatrixSolver( int n,
         double* Diagonal, double* LowerDiag, double* UpperDiag, double* f, double* x )
         {
    double* nDiagonal = Diagonal;
    double* nLowerDiag = LowerDiag;
    double* nUpperDiag = UpperDiag;
    double* nf = f;
                  //scaling matrix entries
    for(int i=1; i<n; i++)
    {
        nDiagonal[i] = nDiagonal[i]-(nUpperDiag[i-1]*nLowerDiag[i]/nDiagonal[i-1]);
        nf[i] = nf[i]-(nf[i-1]*nLowerDiag[i]/nDiagonal[i-1]);
    }

        //solving through elimination
     x[n-1]=nf[n-1]/nDiagonal[n-1];
     for(int i=n-2; i>=0; i=i-1)
     {
         x[i] = (nf[i] - nUpperDiag[i]*x[i+1])/(nDiagonal[i]);
     }
         }
