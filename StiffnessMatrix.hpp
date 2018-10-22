#ifndef STIFFNESSMATRIXHEADERREF
#define STIFFNESSMATRIXHEADERREF
#include "TriDiagMatrix.hpp"
#include <iostream>
#include <vector>
#include <array>

class StiffnessMatrix: public TriDiagMatrix
{
public:

    void BuildStiffnessMatrix ( double a, std::vector<double> SpaceMesh );
};


#endif // MASSMATRIXHEADERREF


