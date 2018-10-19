#ifndef HEATEQUATIONHEADERREF
#define HEATEQUATIONHEADERREF
#include <iostream>

class HeatEquation
{
public:
    void SetInitialMesh();
    void SetVariables();

private:
    void TridiagonalMatrixSolver( int n,
         double* Diagonal, double* LowerDiag, double* UpperDiag, double* f, double* x );

};


#endif // HEATEQUATIONHEADERREF
