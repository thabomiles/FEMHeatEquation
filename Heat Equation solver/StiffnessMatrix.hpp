#ifndef STIFFNESSMATRIXHEADERREF
#define STIFFNESSMATRIXHEADERREF
#include "TriDiagMatrix.hpp"
#include "SpaceMesh.hpp"
#include <iostream>
#include <vector>
#include <array>

class StiffnessMatrix: public TriDiagMatrix
{
public:

    void BuildStiffnessMatrix ( double a, SpaceMesh smesh );
};


#endif // MASSMATRIXHEADERREF


