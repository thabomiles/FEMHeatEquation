#ifndef HEATEQUATIONHEADERREF
#define HEATEQUATIONHEADERREF
#include <iostream>
#include <vector>
#include <array>
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"

class HeatEquation
{
public:

    void SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, const std::string outputFileName );

    void Solve();

    void AnalyticSolutionVec( );

    double PiecewiseU( double x, SpaceMesh currentsmesh, std::vector<double> U  );

    void PrintSolution();

    double ContinuousAnalyticSolution( double x, double t );

    double ErrorSquared( double x );

    double L2ErrorGuass ( double lowerlimit, double upperlimit );

    const double M_PI = 2*acos(0);

    void BuildErrorMesh();

    double GlobalSpaceError();

    void PrintErrorMesh();

    void PrintVector( std::vector<double> aVector);

    void AddVectors( std::vector<double>func1, std::vector<double> func2, std::vector<double>& result );
    void VectorTimesScalar( std::vector<double>& func1, double scalar);

    std::vector<double> mpErrorMesh;



protected:

    int mn, mm, mpcurrenTimeStep = 0, mpcurrentMeshIndex = 0;

    double mT, g_0 =0 , g_L = 0;


    std::string moutputFileName;

    SpaceMesh mpsmesh;
    TimeMesh mptmesh;

    StiffnessMatrix stiff;
    MassMatrix mass;
    TriDiagMatrix LHS;

    std::vector<double> mpx;
    std::vector<double> mpAnalyticSolution;
    std::vector<double> mpPreviousSolution;
    std::vector<double> mpRHS;
    std::vector<double> HalfTimeUh;
    std::vector<double> HalfTimeInterpolant;



//    std::array<double, 2> firstpoint;
//    std::array<double, 2> secondpoint;
};


#endif // HEATEQUATIONHEADERREF
