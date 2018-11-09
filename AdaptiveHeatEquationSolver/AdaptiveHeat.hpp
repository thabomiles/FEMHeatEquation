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
    void SolveChangingMesh();
    void BuildSystemAtTimeStep();
    void SaveIntervalsForRefinement();
    void RefineMesh();
    void UpdatePreviousSolution();
    void SolveStep();
    double InterpolantFunction( double x, std::vector<double> funct, SpaceMesh& relevantMesh );
    void SystemSolver();
    double IntegrateBasisWithU( int NodeIndex, double lowerlimit, double upperlimit, SpaceMesh& currentSmesh,
                               SpaceMesh& previousSmesh, std::vector<double>& SolutionVec );

    double SolutionTimesBasis( int NodeIndex, double x, SpaceMesh& currentSmesh,
                              SpaceMesh& previousSmesh, std::vector<double>& SolutionVec );
    void BuildRHS();




protected:

    const double tolerance = 0.2;
    std::vector<int> intervalsForRefinement;
    SpaceMesh oldmesh;
    SpaceMesh refinedsmesh;
    int mpcurrentMeshIndex = 0;

    std::vector<double> mpPreviousSolDummy;



};


#endif // HEATEQUATIONHEADERREF
