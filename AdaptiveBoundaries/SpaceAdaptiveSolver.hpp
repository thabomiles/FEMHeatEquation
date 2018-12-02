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
    void AdaptiveSolve();
    void SaveIntervalsForRefinement();
    void RefineMesh();

    double IntegrateBasisWithU( int NodeIndex, double lowerlimit, double upperlimit );

    double SolutionTimesBasis( int NodeIndex, double x );
    void BuildRHS();

    void SaveIntervalsForCoarsening();

    void UnitTest();

protected:
    const double tolerance = 0.1, coarseningtol = 0.0;
    std::vector<int> intervalsForRefinement, NodesForRemoval;
    SpaceMesh oldmesh, refinedsmesh;
};


#endif // ADAPTIVESOLVERHEADERREF
