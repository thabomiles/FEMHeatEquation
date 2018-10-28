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
    //number of time steps i.e. 1 less than the number of nodes
int m = mptmesh.NumberOfTimeSteps();

double a = pow(M_PI,-2);

       //n is the number of elements i.e. 1 less than the number of nodes
int n=mpsmesh.meshsize();

n=mpsmesh.meshsize();

StiffnessMatrix stiff;
stiff.BuildStiffnessMatrix( a, mpsmesh );

stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(0) );

MassMatrix mass;
mass.BuildMassMatrix(mpsmesh);
TriDiagMatrix LHS;
LHS.AddTwoMatrices( mass, stiff );


for (int i = 0; i<n-1; i++)
{
    mpPreviousSolution.push_back(6*sin(M_PI*mpsmesh.ReadSpaceNode(i+1)));
}

ofstream myfile;
  myfile.open ("soultion1.txt");

for(int j = 0; j<m; j++)
{
    mpcurrenTimeStep = j+1;

//    mpAnalyticSolution.clear();
//    for (int i = 0; i<n-1; i++)
//{
//    mpAnalyticSolution.push_back(exp(-mptmesh.ReadTimeStep(j+1))*6*
//                               sin(M_PI*mpsmesh.ReadSpaceNode(i+1)));
//}

AnalyticSolution();
mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

LHS.MatrixSolver( mpRHS, mpx );


//        PrintSolution();

        for (auto k: mpx)
        myfile << k <<  " ,"<<' ' ;
        myfile << "\n";


        mpPreviousSolution = mpx;
}
        myfile.close();
}


void HeatEquation::AnalyticSolution( )
{
    mpAnalyticSolution.clear();
     for (int i = 0; i<mpsmesh.meshsize()-1; i++)
{
    mpAnalyticSolution.push_back(exp(-mptmesh.ReadTimeStep(mpcurrenTimeStep))*6*
                               sin(M_PI*mpsmesh.ReadSpaceNode(i+1)));
}
}

void HeatEquation::PrintNodalErrors( )
{
    double nodalerror;
    double infnorm = 0;
    for (int i = 0; i<mpsmesh.meshsize()-1; i++)
{
    nodalerror = fabs(mpAnalyticSolution.at(i)-mpx.at(i));
//    std::cout << nodalerror << ' ';
    if (infnorm<nodalerror)
    {
        infnorm = nodalerror;
    }
}
    std::cout << infnorm << ' ';
    std::cout << " \n";
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



void HeatEquation::PrintSolution( )
{
        for (auto k: mpx)
        std::cout << k << ", ";
        std::cout << " \n";

//        for (auto k: mpAnalyticSolution)
//        std::cout << k << ' ';
//        std::cout << " \n";

}
