#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "heatEquation.hpp"
#include "AdaptiveHeat.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

void AdaptiveHeatEquation::AdaptiveSolver()
{
AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;
BuildSystemAtTimeStep();
oldmesh.CopySpaceMesh(mpsmesh);
mpsmesh.PrintSpaceNodes();


//int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<mptmesh.NumberOfTimeSteps(); j++)
{
    mpcurrenTimeStep = j+1;
    mpcurrentMeshIndex = j;

    BuildSystemAtTimeStep();
    SystemSolver();
    mpPreviousSolution = mpx;

    SaveIntervalsForRefinement();
    SaveIntervalsForCoarsening();
    mpsmesh.BisectIntervals(intervalsForRefinement);
    mpsmesh.CoarsenIntervals(NodesForRemoval);
}
    std::cout<<mpsmesh.meshsize() <<"\n";
    std::cout<<"\n";

}

void AdaptiveHeatEquation::SaveIntervalsForCoarsening()
{
    BuildErrorMesh();
    NodesForRemoval.clear();
    for(int i=0; i<mpErrorMesh.size()-1; i++)
    {
        if (sqrt(mpErrorMesh.at(i)+mpErrorMesh.at(i+1))<coarseningtol)
        {
            NodesForRemoval.push_back(i+1);
        }
    }

}

void AdaptiveHeatEquation::BuildSystemAtTimeStep()
{
stiff.BuildStiffnessMatrix( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );
mass.BuildMassMatrix(mpsmesh);
LHS.AddTwoMatrices( mass, stiff );
}

void AdaptiveHeatEquation::SystemSolver()
{
    BuildRHS();
    LHS.MatrixSolver( mpRHS, mpx );
    oldmesh.CopySpaceMesh(mpsmesh);
}


void AdaptiveHeatEquation::SaveIntervalsForRefinement()
{
    BuildErrorMesh();
    intervalsForRefinement.clear();
    for(int i=0; i<mpErrorMesh.size(); i++)
    {
        if (sqrt(mpErrorMesh.at(i))>tolerance)
        {
            intervalsForRefinement.push_back(i);
        }
    }
}

double AdaptiveHeatEquation::IntegrateBasisWithU( int NodeIndex, double lowerlimit,
                              double upperlimit )
{
    auto SolutionWithBasis = [&](double x)
        { return mpsmesh.TestFunctions( NodeIndex, x)*PiecewiseU(x, oldmesh, mpPreviousSolution); };

    return gauss<double, 7>::integrate(SolutionWithBasis, lowerlimit, upperlimit);
}


double AdaptiveHeatEquation::IntegrateFullBasisWithU( int NodeIndex, double lowerlimit,
                              double upperlimit )
{
   //std::cout<<NodeIndex<<", "<<lowerlimit<<", "<< upperlimit<< "\n";
    auto SolutionWithBasis = [&](double x)
        { return mpsmesh.GeneralTestFunctions( NodeIndex, x)*PiecewiseU(x, oldmesh, mpPreviousSolution); };

    return gauss<double, 7>::integrate(SolutionWithBasis, lowerlimit, upperlimit);
}

    //this function works for refinement and coarsening
void AdaptiveHeatEquation::BuildRHS()
{
mpRHS.clear();
std::vector<double> intervals;

refinedsmesh.CommonMesh(mpsmesh, oldmesh);

double integral;
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
}

void AdaptiveHeatEquation::BuildRHSWithBC()
{
mpRHS.clear();
std::vector<double> intervals;

refinedsmesh.CommonMesh(mpsmesh, oldmesh);

double integral;
refinedsmesh.Range(mpsmesh.ReadSpaceNode(0), mpsmesh.ReadSpaceNode(1), intervals);
    for(int j = 0; j<intervals.size()-1; j++)
    {
        integral = integral + IntegrateFullBasisWithU(0, intervals.at(j), intervals.at(j+1));
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
PrintVector(intervals);
    for(int j = 0; j<intervals.size()-1; j++)
    {
        integral = integral + IntegrateFullBasisWithU(mpsmesh.meshsize(), intervals.at(j), intervals.at(j+1));
    }
    mpRHS.push_back(integral);
}

void AdaptiveHeatEquation::UnitTest()
{
mpsmesh.GenerateUniformMesh(1,4);

oldmesh.GenerateUniformMesh(1, 4);

std::vector<double> u = {2,2};
mpPreviousSolution = u;
std::vector<double> b;
std::vector<double> intervals;

refinedsmesh.CommonMesh(mpsmesh, oldmesh);

MassMatrix test;
test.BuildMassMatrix( mpsmesh );
test.MatrixVectorMultiplier( u ,b);
PrintVector(b);
BuildRHS();
PrintVector(mpRHS);

 u = {0,2,2, 0};
//mpPreviousSolution = u;

test.BuildGeneralMassMatrix(mpsmesh);
test.MatrixVectorMultiplier( u ,b);

PrintVector(b);
b.clear();


BuildRHSWithBC();
PrintVector(mpRHS);
}

