#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "GeneralizedHeat.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

void GeneralHeat::SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, APDE& apde )
{
    mpsmesh = smesh;
    mptmesh = tmesh;
    mppde = &apde;
    k_0 = mppde->k_0;
    k_L = mppde->k_L;
    g_0 = mppde->g_0;
    g_L = mppde->g_L;
}

double GeneralHeat::ContinuousAnalyticSolution( double x, double t )
{
     return mppde->ContinuousAnalyticSolution( x, t );
}

void GeneralHeat::SolveWithBCs()
{
stiff.BuildStiffnessMatrix( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(0) );
mass.BuildMassMatrix(mpsmesh);
LHS.AddTwoMatrices( mass, stiff );

AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
mpcurrenTimeStep = j+1;
stiff.BuildStiffnessMatrix( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrenTimeStep-1) );
mass.BuildMassMatrix(mpsmesh);
LHS.AddTwoMatrices( mass, stiff );

mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

LHS.MatrixSolver( mpRHS, mpx );

mpPreviousSolution = mpx;

if (j==int(0.5*m))
{
    PrintErrorMesh();
    GlobalSpaceError();
}

}
}
