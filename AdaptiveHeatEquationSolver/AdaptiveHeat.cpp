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
    LHS.MatrixSolver( mpRHS, mpx );
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
    PrintVector(mpPreviousSolution);
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

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
    mpcurrenTimeStep = j+1;
    mpcurrentMeshIndex = j;
    BuildSystemAtTimeStep();
    SystemSolver();
    mpsmesh.PrintSpaceNodes();
    PrintSolution();
    mpPreviousSolution = mpx;

if (j==int(0.5*m))
{
    PrintErrorMesh();
    oldmesh.CopySpaceMesh(mpsmesh);
    mpsmesh.GloballyBisectSpaceMesh();
    UpdatePreviousSolution();
    mpsmesh.PrintSpaceNodes();
    std::cout<<"\n";
    std::cout<<mpcurrenTimeStep;
    std::cout<<"\n";
}
}
}

