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
    void BuildSystemAtTimeStep();
    void SaveIntervalsForRefinement();
    void RefineMesh();
    void SolveStep();
    void SystemSolver();
    double IntegrateBasisWithU( int NodeIndex, double lowerlimit, double upperlimit );

    double SolutionTimesBasis( int NodeIndex, double x );
    void BuildRHS();

    void SaveIntervalsForCoarsening();

    void BuildErrorMesh();

    void UnitTest();

    double GlobalSpaceError();


protected:
    const double tolerance = 0.05;
    const double coarseningtol = 0.02;
    std::vector<int> intervalsForRefinement;
    std::vector<int> NodesForRemoval;
    SpaceMesh oldmesh;
    SpaceMesh refinedsmesh;
    int mpcurrentMeshIndex = 0;

};


#endif // ADAPTIVESOLVERHEADERREF
