#ifndef GENERALHEATHEADERREF
#define GENERALHEATHEADERREF
#include <iostream>
#include <vector>
#include <array>
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


protected:
    double k_0 = 0, k_L = 0;
    APDE* mppde;
    std::vector<double> br;

};


#endif // GENERALHEATHEADERREF

