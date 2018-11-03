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
    mpcurrenTimeStep = j+1;

mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

LHS.MatrixSolver( mpRHS, mpx );

//SaveIntervalsForRefinement();

//mpsmesh.InsertSpaceNode(j/m);


mpPreviousSolution = mpx;

//UpdatePreviousSolution();

if (j==int(0.5*m))
{
    PrintSolution();
    BuildErrorMesh();
    PrintErrorMesh();
    std::cout<< GlobalSpaceError();
    std::cout << " \n";
//    mpsmesh.GloballyBisectSpaceMesh();
    UpdatePreviousSolution();
}


}
}

void AdaptiveHeatEquation::UpdatePreviousSolution()
{
//    mpPreviousSolution.clear();
//    double dummyU;
//    double oldMeshSize = mpsmesh.meshsize();
//    for (int i = 1; i<oldMeshSize-2; i++)
//    {
//        dummyU = PiecewiseU(mpsmesh.ReadSpaceNode(i));
//        mpPreviousSolution.push_back(dummyU);
//    }
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

AnalyticSolutionVec();

}

double AdaptiveHeatEquation::InterpolantFunction( double x, std::vector<double> funct )
{
    std::array<double, 2> firstpoint;
    std::array<double, 2> secondpoint;

    int upperindex = mpsmesh.IndexAbove( x );
    auto boundaryconditionU0 = 0;
    auto boundarycondition1Un = 0;

    if((upperindex==1)||(upperindex==0))
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(0);
    firstpoint.at(1) = boundaryconditionU0;

    secondpoint[0] = mpsmesh.ReadSpaceNode(1);
    secondpoint.at(1) = funct.at(0);
    }
    else if (upperindex == mpsmesh.meshsize())
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = funct.at(upperindex-2);

    secondpoint[0] = mpsmesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = boundarycondition1Un;
    }
    else
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = funct.at(upperindex-2);

    secondpoint[0] = mpsmesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = funct.at(upperindex-1);
    }

    long double m = (firstpoint[1]-secondpoint[1])/(firstpoint[0]-secondpoint[0]);

    return m*(x - firstpoint[0])+firstpoint[1];
}



