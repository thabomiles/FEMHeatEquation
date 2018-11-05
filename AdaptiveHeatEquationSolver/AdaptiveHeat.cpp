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
mass.PrintMatrix();
PrintSolution();

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<1; j++)
{
    mpcurrenTimeStep = j+1;
    mpcurrentMeshIndex = j;
    SystemSolver();

    PrintSolution();

    SaveIntervalsForRefinement();
    oldmesh.CopySpaceMesh(mpsmesh);
    mpPreviousSolDummy = mpPreviousSolution;

    int i = 0;
    while (i<1)
    {
        i++;
        mpsmesh.BisectIntervals(intervalsForRefinement);

        UpdatePreviousSolution();
//
        BuildSystemAtTimeStep();
        SystemSolver();
        SaveIntervalsForRefinement();

        if (i==3)
        {
            mpsmesh.PrintSpaceNodes();
            PrintSolution();
        }
    }
    mpsmesh.PrintSpaceNodes();
    PrintSolution();

    mpPreviousSolution = mpx;
}
}

void AdaptiveHeatEquation::UpdatePreviousSolution()
{
//    std::vector<double> previousDummy = mpPreviousSolution;
    mpPreviousSolution.clear();
    double dummyU;
    double dummyx;

    for (int i = 1; i<mpsmesh.meshsize(); i++)
    {
        dummyU = InterpolantFunction(mpsmesh.ReadSpaceNode(i), mpPreviousSolDummy, oldmesh);
        mpPreviousSolution.push_back(dummyU);
    }
//    for (auto i: mpPreviousSolution)
//        std::cout<< i<< ", ";
//    std::cout<< "\n";
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
//    for (auto j: intervalsForRefinement)
//        std::cout<< j << ", ";
//    std::cout<< "\n";
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
BuildSystemAtTimeStep();

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
    mpcurrenTimeStep = j+1;
    mpcurrentMeshIndex = j;

    SystemSolver();

if (j==int(0.5*m))
{
    PrintErrorMesh();
    std::cout << GlobalSpaceError();
    std::cout << " \n";

    oldmesh.CopySpaceMesh(mpsmesh);
    mpsmesh.GloballyBisectSpaceMesh();
    UpdatePreviousSolution();

    BuildSystemAtTimeStep();
    mpsmesh.PrintSpaceNodes();
    PrintSolution();
}
}
}

