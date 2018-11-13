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

double GeneralHeat::ErrorSquared( double x )
{
    double dummyVar = PiecewiseU(x)-ContinuousAnalyticSolution(x, mptmesh.ReadTimeStep(mpcurrenTimeStep));
    return pow(dummyVar,2);
}

double GeneralHeat::PiecewiseU( double x )
{
    std::array<double, 2> firstpoint;
    std::array<double, 2> secondpoint;

    int upperindex = mpsmesh.IndexAbove( x );
    auto boundaryconditionU0 = g_0;
    auto boundarycondition1Un = g_L;

    if((upperindex==1)||(upperindex==0))
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(0);
    firstpoint.at(1) = mpx.at(0);

    secondpoint[0] = mpsmesh.ReadSpaceNode(1);
    secondpoint.at(1) = mpx.at(1);
    }
    else if (upperindex == mpsmesh.meshsize())
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = mpx.at(upperindex-1);

    secondpoint[0] = mpsmesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = mpx.at(upperindex);
    }
    else
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = mpx.at(upperindex-1);

    secondpoint[0] = mpsmesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = mpx.at(upperindex);
    }

    long double m = (firstpoint[1]-secondpoint[1])/(firstpoint[0]-secondpoint[0]);

    return m*(x - firstpoint[0])+firstpoint[1];
}

double GeneralHeat::L2ErrorGuass ( double lowerlimit, double upperlimit )
{
    const int n = 7;
    double halfinterval = (upperlimit-lowerlimit)*pow(2,-1);
    double intervalmidpoint = (upperlimit+lowerlimit)*pow(2,-1);
    auto x  = gauss<double, n>::abscissa();
    auto weight = gauss<double, n>::weights();

    double quad = weight[0]*ErrorSquared(halfinterval*x[0]+intervalmidpoint);

    for (int j = 1; j<=(n-1)*pow(2, -1); j++)
    {
        quad = quad+weight[j]*ErrorSquared(halfinterval*x[j]+intervalmidpoint);
        quad = quad+weight[j]*ErrorSquared(-halfinterval*x[j]+intervalmidpoint);
    }

    return halfinterval*quad;
}

void GeneralHeat::BuildErrorMesh()
{
mpErrorMesh.clear();
for(int i=0; i<mpsmesh.meshsize(); i++)
{
mpErrorMesh.push_back(L2ErrorGuass( mpsmesh.ReadSpaceNode(i),mpsmesh.ReadSpaceNode(i+1)));
}
}

double GeneralHeat::GlobalSpaceError()
{
    BuildErrorMesh();
    double globalError=0;
    for(auto k: mpErrorMesh)
        globalError = globalError + k;

    std::cout << sqrt(globalError);
    std::cout << " \n";
    return sqrt(globalError);
}


void GeneralHeat::PrintErrorMesh()
{
    BuildErrorMesh();
    PrintVector(mpErrorMesh);
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

LHS.AddTwoMatrices( mass, stiff );

mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

BuiltbrVec();
VectorTimesScalar( br, mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );

AddVectors( br, mpRHS, mpRHS );


LHS.MatrixSolver( mpRHS, mpx );

mpPreviousSolution = mpx;


if (j==int(0.5*m))
{
    PrintSolution();
    GlobalSpaceError();
    PrintVector(mpErrorMesh);
    std::cout<<mpcurrenTimeStep;
    std::cout<<"\n";
}

}
BuildErrorMesh();
PrintVector(mpErrorMesh);
}


