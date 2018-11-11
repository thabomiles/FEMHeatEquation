#ifndef MASSMATRIXHEADERREF
#define MASSMATRIXHEADERREF
#include "TriDiagMatrix.hpp"
#include "SpaceMesh.hpp"
#include <iostream>
#include <vector>
#include <array>

class MassMatrix: public TriDiagMatrix
{
public:

    void BuildMassMatrix ( SpaceMesh smesh );
    void BuildGeneralMassMatrix ( SpaceMesh smesh );
};


#endif // MASSMATRIXHEADERREF


