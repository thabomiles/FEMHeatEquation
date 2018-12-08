#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
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

void AdaptiveSolver::SetTolerances(const double refineby, const double coarsenby)
{
    tolerance = refineby;
    coarseningtol = coarsenby;
}

void AdaptiveSolver::AdaptiveChangingBC();
{

}

void AdaptiveSolver::AdaptiveSolve()
{
mpsmesh.PrintSpaceNodes();
mpcurrenTimeStep = 0;
mpcurrentMeshIndex = 0;
oldmesh.CopySpaceMesh(mpsmesh);

AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;
stiff.SetParameters(k_0, k_L, mpa);

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
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
mpsmesh.PrintSpaceNodes();
PrintVector(mpx);
mpPreviousSolution = mpx;

    BuildGradientVec(mpx, mpsmesh, FEMGradient);
    GradientRecoveryFunction( mpsmesh, FEMGradient, GradientRecovery );
    BuildErrorEstimate();

    PrintVector(ErrorEstimate);

    SaveRefinementNodes(mpNodeToInsert);
    SaveIntervalsForCoarsening();
    PrintVector(mpNodeToInsert);

    mpsmesh.CoarsenIntervals(NodesForRemoval);
    mpsmesh.InsertArray( mpNodeToInsert );
}
    std::cout<<mpsmesh.meshsize() <<"\n";
    std::cout<<"\n";
}

void AdaptiveSolver::SaveIntervalsForCoarsening(  )
{
    NodesForRemoval.clear();
    for(int i=0; i<ErrorEstimate.size()-1; i++)
    {
        if (sqrt(ErrorEstimate.at(i)+ErrorEstimate.at(i+1))<coarseningtol)
        {
            NodesForRemoval.push_back(i+1);
        }
    }
}

void AdaptiveSolver::SaveIntervalsForRefinement()
{
    intervalsForRefinement.clear();
    for(int i=0; i<ErrorEstimate.size(); i++)
    {
        if (sqrt(ErrorEstimate.at(i))>tolerance)
        {
            intervalsForRefinement.push_back(i);
        }
    }
}

void AdaptiveSolver::SaveRefinementNodes( std::vector<double>& Nodes_to_insert )
{
    Nodes_to_insert.clear();
    double midpoint;
    for(int i=0; i<ErrorEstimate.size(); i++)
    {
        if (sqrt(ErrorEstimate.at(i))>tolerance)
        {
           midpoint = 0.5*(mpsmesh.ReadSpaceNode(i)+mpsmesh.ReadSpaceNode(i+1));
           Nodes_to_insert.push_back(midpoint);
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

void AdaptiveSolver::MeshRefinement()
{

    SaveIntervalsForRefinement();
    SaveIntervalsForCoarsening();
    mpsmesh.BisectIntervals(intervalsForRefinement);
    mpsmesh.CoarsenIntervals(NodesForRemoval);

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

//void AdaptiveSolver::SaveIntervalsForCoarsening2()
//{
//    NodesForRemoval.clear();
//    for(int i=0; i<ErrorEstimate.size()-1; i++)
//    {
//        if ((sqrt(ErrorEstimate.at(i))<coarseningtol)||(sqrt(ErrorEstimate.at(i+1))<coarseningtol))
//        {
//            NodesForRemoval.push_back(i+1);
//        }
//    }
//}

void AdaptiveSolver::SaveIntervalsForCoarsening2()
{
    NodesForRemoval.clear();
    for(int i=0; i<ErrorEstimate.size()-1; i++)
    {
        if ((sqrt(ErrorEstimate.at(i))<coarseningtol)||(sqrt(ErrorEstimate.at(i+1))<coarseningtol))
        {
            NodesForRemoval.push_back(i+1);
        }
    }
}


