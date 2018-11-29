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

    void AnalyticSolutionVec( );

    void AddVectors( std::vector<double>func1, std::vector<double> func2, std::vector<double>& result );

    void BuildErrorMesh();

    double ContinuousAnalyticSolution( double x, double t );

    double GlobalSpaceError();

    void PrintErrorMesh();

    void PrintVector( std::vector<double> aVector);

    double PiecewiseU( double x, SpaceMesh currentsmesh, std::vector<double> U  );

    void PrintSolution();

    void SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, const std::string outputFileName );

    void Solve();

    void VectorTimesScalar( std::vector<double>& func1, double scalar);

    double GlobalError;
    const double M_PI = 2*acos(0);

protected:

    int mn, mm, mpcurrenTimeStep = 0, mpcurrentMeshIndex = 0;

    double mT;


    std::string moutputFileName;

    SpaceMesh mpsmesh;
    TimeMesh mptmesh;

    StiffnessMatrix stiff;
    MassMatrix mass;
    TriDiagMatrix LHS;

    std::vector<double> mpx, mpAnalyticSolution, mpPreviousSolution, mpRHS, mpErrorMesh;

};


#endif // HEATEQUATIONHEADERREF
