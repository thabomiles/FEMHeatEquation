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
//    std::cout<<"\n";
//    std::cout<<mpcurrenTimeStep <<"\n";
//    std::cout<<"\n";
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

//    for (auto k: NodesForRemoval)
//        std::cout << k << ", ";
//    std::cout << " \n";
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
    //mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );
    BuildRHS();
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
//    for (auto j: mpErrorMesh)
//        std::cout<< sqrt(j) << ", ";
//    std::cout<< "\n";
//
//    for (auto k: intervalsForRefinement)
//        std::cout << k << ", ";
//    std::cout << " \n";
}


double AdaptiveHeatEquation::InterpolantFunction( double x, std::vector<double> funct, SpaceMesh& relevantMesh )
{
    std::array<double, 2> firstpoint;
    std::array<double, 2> secondpoint;

    int upperindex = relevantMesh.IndexAbove( x );
    auto boundaryconditionU0 = g_0;
    auto boundarycondition1Un = g_L;

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

    mpPreviousSolution = mpx;


if (j==int(0.5*m))
{
    //oldmesh.CopySpaceMesh(mpsmesh);
    mpsmesh.GloballyBisectSpaceMesh();
}
}
}

double AdaptiveHeatEquation::IntegrateBasisWithU( int NodeIndex, double lowerlimit,
                              double upperlimit )
{
        const int n = 7;
    double halfinterval = (upperlimit-lowerlimit)*pow(2,-1);
    double intervalmidpoint = (upperlimit+lowerlimit)*pow(2,-1);
    auto x  = gauss<double, n>::abscissa();
    auto weight = gauss<double, n>::weights();

    double quad = weight[0]*SolutionTimesBasis(NodeIndex, halfinterval*x[0]+intervalmidpoint);

    for (int j = 1; j<=(n-1)*pow(2, -1); j++)
    {
        quad = quad+weight[j]*SolutionTimesBasis(NodeIndex, halfinterval*x[j]+intervalmidpoint);
        quad = quad+weight[j]*SolutionTimesBasis(NodeIndex, -halfinterval*x[j]+intervalmidpoint);
    }

    return halfinterval*quad;
}

double AdaptiveHeatEquation::SolutionTimesBasis( int NodeIndex, double x )
{
    return mpsmesh.TestFunctions( NodeIndex, x)*InterpolantFunction(x, mpPreviousSolution, oldmesh);
}

    //this function works for refinement
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




