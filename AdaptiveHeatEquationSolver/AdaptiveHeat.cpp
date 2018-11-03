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

//StiffnessMatrix stiff;
//MassMatrix mass;
//TriDiagMatrix LHS;

//stiff.BuildStiffnessMatrix( mpsmesh );
//stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrenTimeStep) );
//mass.BuildMassMatrix(mpsmesh);
//LHS.AddTwoMatrices( mass, stiff );

BuildSystemAtTimeStep();


AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
    mpcurrenTimeStep = j+1;

AnalyticSolutionVec();
mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

LHS.MatrixSolver( mpRHS, mpx );

mpPreviousSolution = mpx;

if (j==int(0.5*m))
{
    PrintSolution();
    BuildErrorMesh();
    PrintErrorMesh();
    std::cout<< GlobalSpaceError();
    std::cout << " \n";
}

}
}

void AdaptiveHeatEquation::BuildSystemAtTimeStep()
{
stiff.BuildStiffnessMatrix( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrenTimeStep) );
mass.BuildMassMatrix(mpsmesh);
LHS.AddTwoMatrices( mass, stiff );

AnalyticSolutionVec();

}
