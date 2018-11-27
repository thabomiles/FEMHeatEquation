#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "heatEquation.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

void HeatEquation::SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, const std::string outputFileName )
{
    mpsmesh = smesh;
    mptmesh = tmesh;
    moutputFileName = outputFileName;
    mpcurrenTimeStep = 0;
}

void HeatEquation::Solve()
{
AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;

int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
mpcurrenTimeStep = j+1;
mpcurrentMeshIndex = j;
stiff.BuildStiffnessMatrix( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );
mass.BuildMassMatrix(mpsmesh);
LHS.AddTwoMatrices( mass, stiff );

mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

LHS.MatrixSolver( mpRHS, mpx );

mpPreviousSolution = mpx;

if (j==int(0.5*m))
{
//    PrintErrorMesh();
//    GlobalSpaceError();
}

}
}




void HeatEquation::AnalyticSolutionVec( )
{
    mpAnalyticSolution.clear();
     for (int i = 0; i<mpsmesh.meshsize()-1; i++)
{

    mpAnalyticSolution.push_back(ContinuousAnalyticSolution( mpsmesh.ReadSpaceNode(i+1),
                                                                        mptmesh.ReadTimeStep(mpcurrenTimeStep)));
}
}


    //the if statements deal with what happens if you are in the first interval
    //in this case the first value of u is given by the boundary condition
    //all the other elements can be fixed for x = 0
    //this function only works for the static case if the mesh changes everything is fucked
double HeatEquation::PiecewiseU( double x, SpaceMesh currentsmesh, std::vector<double> U )
{
    std::array<double, 2> firstpoint;
    std::array<double, 2> secondpoint;

    int upperindex = currentsmesh.IndexAbove( x );
    auto boundaryconditionU0 = g_0;
    auto boundarycondition1Un = g_L;

    if((upperindex==1)||(upperindex==0))
    {
    firstpoint.at(0)= currentsmesh.ReadSpaceNode(0);
    firstpoint.at(1) = boundaryconditionU0;

    secondpoint[0] = currentsmesh.ReadSpaceNode(1);
    secondpoint.at(1) = U.at(0);
    }
    else if (upperindex == currentsmesh.meshsize())
    {
    firstpoint.at(0)= currentsmesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = U.at(upperindex-2);

    secondpoint[0] = currentsmesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = boundarycondition1Un;
    }
    else
    {
    firstpoint.at(0)= currentsmesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = U.at(upperindex-2);

    secondpoint[0] = currentsmesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = U.at(upperindex-1);
    }

    long double m = (firstpoint[1]-secondpoint[1])/(firstpoint[0]-secondpoint[0]);

    return m*(x - firstpoint[0])+firstpoint[1];
}


double HeatEquation::ContinuousAnalyticSolution( double x, double t )
{
     return 6*exp(-t)*sin(M_PI*x);
}


void HeatEquation::BuildErrorMesh()
{
mpErrorMesh.clear();

auto SquaredError = [&](double x)
    { return pow(PiecewiseU(x, mpsmesh, mpx)-ContinuousAnalyticSolution(x, mptmesh.ReadTimeStep(mpcurrenTimeStep)), 2); };

double Q;
for(int i=0; i<mpsmesh.meshsize(); i++)
{
Q = gauss<double, 7>::integrate(SquaredError, mpsmesh.ReadSpaceNode(i), mpsmesh.ReadSpaceNode(i+1));
mpErrorMesh.push_back( Q );
}
}

double HeatEquation::GlobalSpaceError()
{
    BuildErrorMesh();
    double globalError=0;
    for(auto k: mpErrorMesh)
        globalError = globalError + k;

    std::cout << sqrt(globalError);
    std::cout << " \n";
    return sqrt(globalError);
}

void HeatEquation::PrintErrorMesh()
{
    BuildErrorMesh();
    PrintVector(mpErrorMesh);
}


void HeatEquation::PrintSolution( )
{
        AnalyticSolutionVec();
        PrintVector(mpx);
        PrintVector(mpAnalyticSolution);
}

void HeatEquation::PrintVector( std::vector<double> aVector)
{
        for (auto k: aVector)
        std::cout << k << ", ";
        std::cout << " \n";

}

void HeatEquation::AddVectors(std::vector<double>func1, std::vector<double> func2, std::vector<double>& result)
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

void HeatEquation::VectorTimesScalar( std::vector<double>& func1, double scalar)
{
         for ( int i=0; i<func1.size(); i++ )
     {
         func1.at(i)= scalar*func1.at(i);
     }
}



