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
#include <string>
#include <functional>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

void AdaptiveSolver::AdaptiveSolve()
{
oldmesh.CopySpaceMesh(mpsmesh);
mpsmesh.PrintSpaceNodes();

AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;
stiff.SetParameters(k_0, k_L, mpa);


int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
//std::cout<< mpcurrenTimeStep<< ", ";
mpcurrenTimeStep = j+1;
mpcurrentMeshIndex = j;
stiff.BuildGeneralStiffnessMatrix ( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );
mass.BuildGeneralMassMatrix(mpsmesh);

LHS.AddTwoMatrices( mass, stiff );

BuildRHS();
BuiltbrVec();
VectorTimesScalar( br, mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );

AddVectors( br, mpRHS, mpRHS );

LHS.MatrixSolver( mpRHS, mpx );
oldmesh.CopySpaceMesh(mpsmesh);

mpPreviousSolution = mpx;


    SaveIntervalsForRefinement();
    SaveIntervalsForCoarsening();
    mpsmesh.BisectIntervals(intervalsForRefinement);
    mpsmesh.CoarsenIntervals(NodesForRemoval);
}
    std::cout<<mpsmesh.meshsize() <<"\n";
    std::cout<<"\n";
}

void AdaptiveSolver::SaveIntervalsForCoarsening()
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

void AdaptiveSolver::BuildSystemAtTimeStep()
{
stiff.BuildGeneralStiffnessMatrix( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );
mass.BuildGeneralMassMatrix(mpsmesh);
LHS.AddTwoMatrices( mass, stiff );
}

void AdaptiveSolver::SystemSolver()
{
    BuildRHS();
    LHS.MatrixSolver( mpRHS, mpx );
    oldmesh.CopySpaceMesh(mpsmesh);
}


void AdaptiveSolver::SaveIntervalsForRefinement()
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
    for(int j = 0; j<intervals.size()-1; j++)
    {
        integral = integral + IntegrateBasisWithU(mpsmesh.meshsize(), intervals.at(j), intervals.at(j+1));
    }
    mpRHS.push_back(integral);

//PrintVector(mpRHS);
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
mpPreviousSolution = {1,2,2,1};


BuildRHS();

MassMatrix test;
test.BuildGeneralMassMatrix(mpsmesh);
test.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );
PrintVector(mpRHS);

mpsmesh.RemoveSpaceNode(1);
mpsmesh.PrintSpaceNodes();

BuildRHS();
}


void AdaptiveSolver::BuildErrorMesh()
{
mpErrorMesh.clear();

auto SquaredError = [this](double x)
    { return pow(GeneralInterpolant(x, mpx, mpsmesh ) -
    ContinuousAnalyticSolution(x, mptmesh.ReadTimeStep(mpcurrenTimeStep)), 2); };

double Q;
for(int i=0; i<mpsmesh.meshsize(); i++)
{
Q = gauss<double, 7>::integrate(SquaredError, mpsmesh.ReadSpaceNode(i), mpsmesh.ReadSpaceNode(i+1));
mpErrorMesh.push_back( Q );
}
}

double AdaptiveSolver::GlobalSpaceError()
{
    BuildErrorMesh();
    double globalError=0;
    for(auto k: mpErrorMesh)
        globalError = globalError + k;

    std::cout << sqrt(globalError);
    std::cout << " \n";
    return sqrt(globalError);
}
