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
    void Solve();
    void BuildSystemAtTimeStep();
    void SaveIntervalsForRefinement();
    void RefineMesh();
    void UpdatePreviousSolution();
    void SolveStep();
    double InterpolantFunction( double x, std::vector<double> funct, SpaceMesh& relevantMesh );
    void SolveTimeStep();




protected:

    const double tolerance = 0.1;
    std::vector<double> intervalsForRefinement;
    SpaceMesh oldmesh;



};


#endif // HEATEQUATIONHEADERREF
