#ifndef STIFFNESSMATRIXHEADERREF
#define STIFFNESSMATRIXHEADERREF
#include "TriDiagMatrix.hpp"
#include "SpaceMesh.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <math.h>

class StiffnessMatrix: public TriDiagMatrix
{
public:

    void BuildStiffnessMatrix ( SpaceMesh smesh );

protected:
    const double M_PI = 2*acos(0);
    double a = pow(M_PI,-2);
};


#endif // MASSMATRIXHEADERREF


