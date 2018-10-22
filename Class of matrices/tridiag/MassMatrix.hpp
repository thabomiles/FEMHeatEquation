#ifndef MASSMATRIXHEADERREF
#define MASSMATRIXHEADERREF
#include "TriDiagMatrix.hpp"
#include <iostream>
#include <vector>
#include <array>

class MassMatrix: public TriDiagMatrix
{
public:

    void BuildMassMatrix ( std::vector<double> SpaceMesh );
};


#endif // MASSMATRIXHEADERREF


