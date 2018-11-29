#ifndef ADAPTIVEHEATEQUATIONHEADERREF
#define ADAPTIVEHEATEQUATIONHEADERREF
#include <iostream>
#include <vector>
#include <array>
#include "heatEquation.hpp"
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"

class AdaptiveHeatEquation: public HeatEquation
{
public:
    void AdaptiveSolver();
    void BuildSystemAtTimeStep();
    void SaveIntervalsForRefinement();
    void RefineMesh();
    void SolveStep();
    void SystemSolver();
    double IntegrateBasisWithU( int NodeIndex, double lowerlimit, double upperlimit );

    double SolutionTimesBasis( int NodeIndex, double x );
    void BuildRHS();

    void SaveIntervalsForCoarsening();

    void UnitTest();



protected:

    const double tolerance = 0.2, coarseningtol = 0.05;
    std::vector<int> intervalsForRefinement, NodesForRemoval;
    SpaceMesh oldmesh, refinedsmesh;
    int mpcurrentMeshIndex = 0;

};


#endif // ADAPTIVEHEATEQUATIONHEADERREF
