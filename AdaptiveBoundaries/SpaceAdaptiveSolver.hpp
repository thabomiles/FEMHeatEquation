#ifndef ADAPTIVESOLVERHEADERREF
#define ADAPTIVESOLVERHEADERREF
#include <iostream>
#include <vector>
#include <array>
#include <functional>
#include "GeneralizedHeat.hpp"
#include "AdaptiveHeat.hpp"
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include "APDE.hpp"

class AdaptiveSolver: public GeneralHeat
{
public:
    double SolutionTimesBasis( int NodeIndex, double x );
    void BuildRHS();
    double IntegrateBasisWithU( int NodeIndex, double lowerlimit, double upperlimit );
    void UnitTest();


protected:

};


#endif // ADAPTIVESOLVERHEADERREF
