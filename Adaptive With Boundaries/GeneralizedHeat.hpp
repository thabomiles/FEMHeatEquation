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

class GeneralHeat: public AdaptiveHeatEquation
{
public:
    void SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, APDE& apde );
    double ContinuousAnalyticSolution( double x, double t );
    void BuiltbrVec();

    void SolveWithBCs();
    void AnalyticSolutionVec( );
    void PrintSolution( );
    double PiecewiseU( double x );

    double ErrorSquared( double x );

    double L2ErrorGuass ( double lowerlimit, double upperlimit );

    void BuildErrorMesh();

    double GlobalSpaceError();

    void PrintErrorMesh();

    double GeneralInterpolant( double x, std::vector<double>& funct, SpaceMesh& relevantMesh );

    void BuildGradientVec( std::vector<double>& funct, SpaceMesh& relevantMesh, std::vector<double>& gradvec );

    void GradientRecoveryFunction(  SpaceMesh& relevantMesh,
                                             std::vector<double>& gradvec, std::vector<double>& gradrecovery );

    void BuildErrorEstimate ( );

    double GradientFunction ( double x );

    void StationaryHeatEquation();

    void UnitTest1 ();

    void buildfvec( SpaceMesh& a_smesh, const std::function<double(double)>& f );

    std::vector<double> FEMGradient;
    std::vector<double> ErrorEstimate;
    std::vector<double> GradientRecovery;

protected:
    double k_0, k_L, g_0, g_L;
    APDE* mppde;
    std::vector<double> br;
    std::vector<double> f_vec;


};


#endif // GENERALHEATHEADERREF

