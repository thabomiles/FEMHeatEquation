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

void AdaptiveHeatEquation::Solve()
{
AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;



int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{

    BuildSystemAtTimeStep();
    oldmesh.CopySpaceMesh(mpsmesh);

    mpcurrenTimeStep = j+1;

    mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );
    LHS.MatrixSolver( mpRHS, mpx );
    mpPreviousSolution = mpx;




if (j==int(0.5*m))
{


    BuildErrorMesh();
    PrintErrorMesh();

    std::cout << GlobalSpaceError();
    std::cout << " \n";

    mpsmesh.GloballyBisectSpaceMesh();
    UpdatePreviousSolution();
    mpsmesh.PrintSpaceNodes();
    PrintSolution();
}
}
}

void AdaptiveHeatEquation::UpdatePreviousSolution()
{
    mpPreviousSolution.clear();
    double dummyU;
    double dummyx;

    for (int i = 1; i<mpsmesh.meshsize(); i++)
    {
        dummyU = InterpolantFunction(mpsmesh.ReadSpaceNode(i), mpx, oldmesh);
        mpPreviousSolution.push_back(dummyU);
    }
//    for (auto i: mpPreviousSolution)
//        std::cout<< i<< ", ";
//    std::cout<< "\n";
}

void AdaptiveHeatEquation::RefineMesh()
{
    for(int i = 0; i<intervalsForRefinement.size(); i++)
    {
        mpsmesh.BisectInterval(i, i+1);
    }
    BuildErrorMesh();
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

    for (auto j: intervalsForRefinement)
        std::cout<< j << ", ";
    std::cout<< "\n";
}

void AdaptiveHeatEquation::BuildSystemAtTimeStep()
{
stiff.BuildStiffnessMatrix( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrenTimeStep) );
mass.BuildMassMatrix(mpsmesh);
LHS.AddTwoMatrices( mass, stiff );
}

void AdaptiveHeatEquation::SolveTimeStep()
{

mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

LHS.MatrixSolver( mpRHS, mpx );

mpPreviousSolution = mpx;

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



