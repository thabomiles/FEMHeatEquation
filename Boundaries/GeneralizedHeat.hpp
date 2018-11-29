#ifndef GENERALHEATHEADERREF
#define GENERALHEATHEADERREF
#include <iostream>
#include <vector>
#include <array>
#include <functional>
#include "AdaptiveHeat.hpp"
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include "APDE.hpp"

class GeneralHeat
{
public:
    void AnalyticSolutionVec( );

    void AddVectors( std::vector<double>func1, std::vector<double> func2, std::vector<double>& result );

    void buildfvec( SpaceMesh& a_smesh );

    void BuiltbrVec();

    void BuildErrorMesh();

    void BuildGradientVec( std::vector<double>& funct, SpaceMesh& relevantMesh, std::vector<double>& gradvec );

    void BuildErrorEstimate ( );

    double ContinuousAnalyticSolution( double x, double t );

    double GlobalSpaceError();

    double GeneralInterpolant( double x, std::vector<double>& funct, SpaceMesh& relevantMesh );

    void GradientRecoveryFunction(  SpaceMesh& relevantMesh,
                                             std::vector<double>& gradvec, std::vector<double>& gradrecovery );

    double GradientFunction ( double x );

    void PrintErrorMesh();

    void PrintSolution( );

    void PrintVector( std::vector<double> aVector);

    void SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, APDE& apde );

    void SolveWithBCs();

    void StationaryHeatEquation();

    void UnitTest1 ();

    void VectorTimesScalar( std::vector<double>& func1, double scalar);


    const double M_PI = 2*acos(0);

    std::vector<double> FEMGradient, ErrorEstimate, GradientRecovery;

protected:
    double mpa,mT, k_0, k_L, g_0, g_L;

    int mn, mm, mpcurrenTimeStep = 0, mpcurrentMeshIndex = 0;

    APDE* mppde;
    std::vector<double> br, f_vec, mpx, mpAnalyticSolution, mpPreviousSolution, mpRHS, mpErrorMesh;

    SpaceMesh mpsmesh;
    TimeMesh mptmesh;

    StiffnessMatrix stiff;
    MassMatrix mass;
    TriDiagMatrix LHS;


};


#endif // GENERALHEATHEADERREF

