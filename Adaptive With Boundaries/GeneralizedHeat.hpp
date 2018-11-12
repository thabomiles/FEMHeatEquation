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
    void AddVectors( std::vector<double>func1, std::vector<double> func2, std::vector<double>& result );
    void VectorTimesScalar( std::vector<double>& func1, double scalar);
    void SolveWithBCs();
    void AnalyticSolutionVec( );
    void PrintSolution( );


protected:
    double k_0 = 0, k_L = 0;
    APDE* mppde;
    std::vector<double> br;

};


#endif // GENERALHEATHEADERREF

