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

    void SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, const std::string outputFileName );

    void Solve();

    void AnalyticSolutionVec( );

    double PiecewiseU( double x );

    void PrintSolution();

    double ContinuousAnalyticSolution( double x, double t );

    double ErrorSquared( double x );

    double L2ErrorGuass ( double lowerlimit, double upperlimit );

    const double M_PI = 2*acos(0);

    void BuildErrorMesh();

    double GlobalSpaceError();

    void PrintErrorMesh();

    std::vector<double> mpErrorMesh;

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


    std::array<double, 2> firstpoint;
    std::array<double, 2> secondpoint;
};


#endif // HEATEQUATIONHEADERREF
