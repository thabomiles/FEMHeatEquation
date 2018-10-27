#ifndef HEATEQUATIONHEADERREF
#define HEATEQUATIONHEADERREF
#include <iostream>
#include <vector>
#include <array>
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"

class HeatEquation
{
public:

    double M_PI = 2*acos(0);

    void SolveUniformMesh();

    void SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, const std::string outputFileName );

    void Solve();

    void AnalyticSolution( );

    void PrintNodalErrors( );



private:

    int mn, mm, mpcurrenTimeStep;

    double mT;

    std::string moutputFileName;

    SpaceMesh mpsmesh;
    TimeMesh mptmesh;

    std::vector<double> mpx;
    std::vector<double> mpAnalyticSolution;
    std::vector<double> mpPreviousSolution;
    std::vector<double> mpRHS;
};


#endif // HEATEQUATIONHEADERREF
