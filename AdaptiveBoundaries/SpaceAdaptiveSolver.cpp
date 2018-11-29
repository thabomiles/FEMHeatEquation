#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "SpaceAdaptiveSolver.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <fstream>
#include <string>
#include <functional>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;


void AdaptiveSolver::BuildRHS()
{
mpRHS.clear();
std::vector<double> intervals;

refinedsmesh.CommonMesh(mpsmesh, oldmesh);

double integral;
refinedsmesh.Range(mpsmesh.ReadSpaceNode(0), mpsmesh.ReadSpaceNode(1), intervals);
    for(int j = 0; j<intervals.size()-1; j++)
    {
        integral = integral + IntegrateBasisWithU(0, intervals.at(j), intervals.at(j+1));
    }
    mpRHS.push_back(integral);

for (int i = 1; i<mpsmesh.meshsize(); i++)
{
    integral=0;
    refinedsmesh.Range(mpsmesh.ReadSpaceNode(i-1), mpsmesh.ReadSpaceNode(i+1), intervals);

    for(int j = 0; j<intervals.size()-1; j++)
    {
        integral = integral + IntegrateBasisWithU(i, intervals.at(j),
                    intervals.at(j+1));
    }
    mpRHS.push_back(integral);
}
integral=0;
refinedsmesh.Range(mpsmesh.ReadSpaceNode(mpsmesh.meshsize()-1), mpsmesh.ReadSpaceNode(mpsmesh.meshsize()), intervals);
//PrintVector(intervals);
    for(int j = 0; j<intervals.size()-1; j++)
    {
        integral = integral + IntegrateBasisWithU(mpsmesh.meshsize(), intervals.at(j), intervals.at(j+1));
    }
    mpRHS.push_back(integral);
}


double AdaptiveSolver::IntegrateBasisWithU( int NodeIndex, double lowerlimit,
                              double upperlimit )
{
   //std::cout<<NodeIndex<<", "<<lowerlimit<<", "<< upperlimit<< "\n";
    auto SolutionWithBasis = [&](double x)
        { return mpsmesh.GeneralTestFunctions( NodeIndex, x)*GeneralInterpolant( x, mpPreviousSolution, oldmesh ); };

    return gauss<double, 7>::integrate(SolutionWithBasis, lowerlimit, upperlimit);
}

void AdaptiveSolver::UnitTest()
{
mpsmesh.GenerateUniformMesh(1,4);

oldmesh.GenerateUniformMesh(1, 4);

std::vector<double> u = {1,2,2,1};
mpPreviousSolution = u;
std::vector<double> b;
std::vector<double> intervals;

refinedsmesh.CommonMesh(mpsmesh, oldmesh);

MassMatrix test;
BuildRHS();
PrintVector(mpRHS);


test.BuildGeneralMassMatrix(mpsmesh);
test.MatrixVectorMultiplier( u ,b);

PrintVector(b);

}


