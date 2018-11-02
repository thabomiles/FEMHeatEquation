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
StiffnessMatrix stiff;
stiff.BuildStiffnessMatrix( mpsmesh );

stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(0) );

MassMatrix mass;
mass.BuildMassMatrix(mpsmesh);
TriDiagMatrix LHS;
LHS.AddTwoMatrices( mass, stiff );


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

}
}


void HeatEquation::AnalyticSolutionVec( )
{
    mpAnalyticSolution.clear();
     for (int i = 0; i<mpsmesh.meshsize()-1; i++)
{

    mpAnalyticSolution.push_back(HeatEquation::ContinuousAnalyticSolution( mpsmesh.ReadSpaceNode(i+1),
                                                                        mptmesh.ReadTimeStep(mpcurrenTimeStep)));
}
}


    //the if statements deal with what happens if you are in the first interval
    //in this case the first value of u is given by the boundary condition
    //all the other elements can be fixed for x = 0
double HeatEquation::PiecewiseU( double x )
{
    int upperindex = mpsmesh.IndexAbove( x );
    auto boundaryconditionU0 = 0;
    auto boundarycondition1Un = 0;

    if((upperindex==1)||(upperindex==0))
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(0);
    firstpoint.at(1) = boundaryconditionU0;

    secondpoint[0] = mpsmesh.ReadSpaceNode(1);
    secondpoint.at(1) = mpx.at(0);
    }
    else if (upperindex == mpsmesh.meshsize())
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = mpx.at(upperindex-2);

    secondpoint[0] = mpsmesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = boundarycondition1Un;
    }
    else
    {
    firstpoint.at(0)= mpsmesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = mpx.at(upperindex-2);

    secondpoint[0] = mpsmesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = mpx.at(upperindex-1);
    }



    double m = (firstpoint[1]-secondpoint[1])/(firstpoint[0]-secondpoint[0]);

    return m*(x - firstpoint[0])+firstpoint[1];

}

double HeatEquation::ErrorSquared( double x )
{
    double dummyVar = HeatEquation::PiecewiseU(x)-HeatEquation::ContinuousAnalyticSolution(x, 1);
    return pow(dummyVar,2);
}

double HeatEquation::ContinuousAnalyticSolution( double x, double t )
{
     return 6*exp(-t)*sin(M_PI*x);
}


double HeatEquation::L2ErrorGuass ( double lowerlimit, double upperlimit )
{
    const int n = 7;
    double halfinterval = (upperlimit-lowerlimit)*pow(2,-1);
    double intervalmidpoint = (upperlimit+lowerlimit)*pow(2,-1);
    auto x  = gauss<double, n>::abscissa();
    auto weight = gauss<double, n>::weights();

    double quad = weight[0]*HeatEquation::ErrorSquared(halfinterval*x[0]+intervalmidpoint);
//    std::cout << quad << ", ";

    for (int j = 1; j<=(n-1)*pow(2, -1); j++)
    {
        quad = quad+weight[j]*HeatEquation::ErrorSquared(halfinterval*x[j]+intervalmidpoint);
        quad = quad+weight[j]*HeatEquation::ErrorSquared(-halfinterval*x[j]+intervalmidpoint);
    }

    return halfinterval*quad;
}

void HeatEquation::BuildErrorMesh()
{
mpErrorMesh.clear();
for(int i=0; i<mpsmesh.meshsize(); i++)
{
mpErrorMesh.push_back(HeatEquation::L2ErrorGuass( mpsmesh.ReadSpaceNode(i),mpsmesh.ReadSpaceNode(i+1)));
}
}

double HeatEquation::GlobalSpaceError()
{
    double globalError=0;
    for(auto k: mpErrorMesh)
        globalError = globalError + k;

    return sqrt(globalError);
}

void HeatEquation::PrintErrorMesh()
{
    for(auto k: mpErrorMesh)
        std::cout << k << ", ";
        std::cout << " \n";

}


void HeatEquation::PrintSolution( )
{
        for (auto k: mpx)
        std::cout << k << ", ";
        std::cout << " \n";

        for (auto k: mpAnalyticSolution)
        std::cout << k << ", ";
        std::cout << " \n";

}
