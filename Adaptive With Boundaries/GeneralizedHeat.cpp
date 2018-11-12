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

void GeneralHeat::BuiltbrVec()
{
    br.assign(mpsmesh.meshsize()+1, 0);
    br.at(0) = k_0*g_0;
    br.at(mpsmesh.meshsize()) = k_L*g_L;


}

void GeneralHeat::AddVectors(std::vector<double>func1, std::vector<double> func2, std::vector<double>& result)
{
    if(func1.size()==func2.size())
    {
    result.clear();
     for ( int i=0; i<func1.size(); i++ )
     {
         result.push_back(func1.at(i)+func2.at(i));
     }
    }
    else
    {
        std::cout<< " your vectors are different sizes";
        std::cout<<"\n";
    }
}

void GeneralHeat::VectorTimesScalar( std::vector<double>& func1, double scalar)
{
         for ( int i=0; i<func1.size(); i++ )
     {
         func1.at(i)= scalar*func1.at(i);
     }
}

void GeneralHeat::AnalyticSolutionVec( )
{
    mpAnalyticSolution.clear();
     for (int i = 0; i<mpsmesh.meshsize()+1; i++)
{

    mpAnalyticSolution.push_back(ContinuousAnalyticSolution( mpsmesh.ReadSpaceNode(i),
                                                                        mptmesh.ReadTimeStep(mpcurrenTimeStep)));
}
}

void GeneralHeat::PrintSolution( )
{
        AnalyticSolutionVec();
        PrintVector(mpx);
        PrintVector(mpAnalyticSolution);
}

void GeneralHeat::SolveWithBCs()
{
AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;
stiff.SetParameters(k_0, k_L);

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
mpcurrenTimeStep = j+1;
mpcurrentMeshIndex = j;
stiff.BuildGeneralStiffnessMatrix ( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );
mass.BuildGeneralMassMatrix(mpsmesh);
//stiff.PrintMatrix();
//std::cout<<"\n";
//mass.PrintMatrix();
LHS.AddTwoMatrices( mass, stiff );
//std::cout<<"\n";
//LHS.PrintMatrix();
//std::cout<<"\n";
mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

BuiltbrVec();
VectorTimesScalar( br, mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );

AddVectors( br, mpRHS, mpRHS );


LHS.MatrixSolver( mpRHS, mpx );

PrintSolution();
std::cout<<"\n";
GlobalSpaceError();

mpPreviousSolution = mpx;

if (j==int(0.5*m))
{
    PrintErrorMesh();
    GlobalSpaceError();
}

}
}
