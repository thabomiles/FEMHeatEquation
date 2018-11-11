#ifndef GENERALHEATHEADERREF
#define GENERALHEATHEADERREF
#include <iostream>
#include <vector>
#include <array>
#include "AdaptiveHeat.hpp"
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include "APDE.hpp"

class GeneralHeat: public AdaptiveHeatEquation
{
public:
    void SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, APDE& apde );
    double ContinuousAnalyticSolution( double x, double t );
    void SolveWithBCs();

protected:
    APDE* mppde;


};


#endif // GENERALHEATHEADERREF

