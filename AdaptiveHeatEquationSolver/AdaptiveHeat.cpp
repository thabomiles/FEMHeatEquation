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

//int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<mptmesh.NumberOfTimeSteps(); j++)
{
    mpcurrenTimeStep = j+1;
    mpcurrentMeshIndex = j;
    BuildSystemAtTimeStep();
    SystemSolver();
    PrintSolution();
    mpPreviousSolution = mpx;

    SaveIntervalsForRefinement();
    oldmesh.CopySpaceMesh(mpsmesh);
    mpsmesh.BisectIntervals(intervalsForRefinement);
    mpsmesh.PrintSpaceNodes();
    UpdatePreviousSolution();
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
    mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );
    //BuildRHS();
    LHS.MatrixSolver( mpRHS, mpx );
    oldmesh.CopySpaceMesh(mpsmesh);
}

void AdaptiveHeatEquation::UpdatePreviousSolution()
{
    mpPreviousSolution.clear();
    double dummyU;

    for (int i = 1; i<mpsmesh.meshsize(); i++)
    {
        dummyU = InterpolantFunction(mpsmesh.ReadSpaceNode(i), mpx, oldmesh);
        mpPreviousSolution.push_back(dummyU);
    }
    //PrintVector(mpPreviousSolution);
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
    for (auto j: mpErrorMesh)
        std::cout<< sqrt(j) << ", ";
    std::cout<< "\n";

    for (auto k: intervalsForRefinement)
        std::cout << k << ", ";
    std::cout << " \n";
}


double AdaptiveHeatEquation::InterpolantFunction( double x, std::vector<double> funct, SpaceMesh& relevantMesh )
{
    std::array<double, 2> firstpoint;
    std::array<double, 2> secondpoint;

    int upperindex = relevantMesh.IndexAbove( x );
    auto boundaryconditionU0 = 0;
    auto boundarycondition1Un = 0;

    if((upperindex==1)||(upperindex==0))
    {
    firstpoint.at(0)= relevantMesh.ReadSpaceNode(0);
    firstpoint.at(1) = boundaryconditionU0;

    secondpoint[0] = relevantMesh.ReadSpaceNode(1);
    secondpoint.at(1) = funct.at(0);
    }
    else if (upperindex == relevantMesh.meshsize())
    {
    firstpoint.at(0)= relevantMesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = funct.at(upperindex-2);

    secondpoint[0] = relevantMesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = boundarycondition1Un;
    }
    else
    {
    firstpoint.at(0)= relevantMesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = funct.at(upperindex-2);

    secondpoint[0] = relevantMesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = funct.at(upperindex-1);
    }

    long double m = (firstpoint[1]-secondpoint[1])/(firstpoint[0]-secondpoint[0]);

    return m*(x - firstpoint[0])+firstpoint[1];
}

void AdaptiveHeatEquation::SolveChangingMesh()
{
AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;
oldmesh.CopySpaceMesh(mpsmesh);

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
    mpcurrenTimeStep = j+1;
    mpcurrentMeshIndex = j;

    BuildSystemAtTimeStep();
    //mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );
    BuildRHS();
    LHS.MatrixSolver( mpRHS, mpx );
    oldmesh.CopySpaceMesh(mpsmesh);
    PrintVector(mpRHS);

    mpPreviousSolution = mpx;


if (j==int(0.5*m))
{
    //oldmesh.CopySpaceMesh(mpsmesh);
    mpsmesh.GloballyBisectSpaceMesh();
}
}
}

double AdaptiveHeatEquation::IntegrateBasisWithU( int NodeIndex, double lowerlimit,
                              double upperlimit, SpaceMesh& currentSmesh, SpaceMesh& previousSmesh,
                               std::vector<double>& SolutionVec )
{
        const int n = 7;
    double halfinterval = (upperlimit-lowerlimit)*pow(2,-1);
    double intervalmidpoint = (upperlimit+lowerlimit)*pow(2,-1);
    auto x  = gauss<double, n>::abscissa();
    auto weight = gauss<double, n>::weights();

    double quad = weight[0]*SolutionTimesBasis(NodeIndex, halfinterval*x[0]+intervalmidpoint,
                                            currentSmesh, previousSmesh, SolutionVec);

    for (int j = 1; j<=(n-1)*pow(2, -1); j++)
    {
        quad = quad+weight[j]*SolutionTimesBasis(NodeIndex, halfinterval*x[j]+intervalmidpoint,
                                                 currentSmesh, previousSmesh, SolutionVec);
        quad = quad+weight[j]*SolutionTimesBasis(NodeIndex, -halfinterval*x[j]+intervalmidpoint,
                                                 currentSmesh, previousSmesh, SolutionVec);
    }

    return halfinterval*quad;
}

double AdaptiveHeatEquation::SolutionTimesBasis( int NodeIndex, double x, SpaceMesh& currentSmesh,
                              SpaceMesh& previousSmesh, std::vector<double>& SolutionVec )
{
    return currentSmesh.TestFunctions( NodeIndex, x)*InterpolantFunction(x, SolutionVec, previousSmesh);
}

    //this function is untested in context
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
                    intervals.at(j+1), mpsmesh, oldmesh, mpPreviousSolution);
    }
    mpRHS.push_back(integral);
}
}




