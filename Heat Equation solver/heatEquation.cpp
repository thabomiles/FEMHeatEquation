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
#include "GeneralMesh.hpp"
#include <fstream>
#include <string>
using namespace std;

void HeatEquation::SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, const std::string outputFileName )
{
    mpSpaceMesh = smesh;
    mpTimeMesh = tmesh;
    moutputFileName = outputFileName;

}

void HeatEquation::Solve()
{

double M_PI = 2*acos(0);

double a = pow(M_PI,-2);


StiffnessMatrix stiff;
stiff.BuildStiffnessMatrix( a, mpSpaceMesh );
stiff.MultiplyByScalar( mpTimeMesh.ReadTimeMesh(0) );


MassMatrix mass;
mass.BuildMassMatrix( mpSpaceMesh );

TriDiagMatrix LHS;
LHS.AddTwoMatrices( mass, stiff );


std::vector<double> x;
std::vector<double> AnalyticSolution;
std::vector<double> PreviousSolution;
std::vector<double> RHS;

for (int i = 0; i<mn-1; i++)
{
    PreviousSolution.push_back(6*sin(M_PI*mpSpaceMesh.ReadSpaceNode(i+1)));
    //PreviousSolution.push_back(6*sin(M_PI*initialSpaceNodes.at(i+1)));
}

ofstream myfile;
  myfile.open (moutputFileName);

for(int j = 0; j<mm; j++)
{
    for (int i = 0; i<mn-1; i++)
{
    AnalyticSolution.push_back(exp(-mpTimeMesh.ReadTimeStep(j+1))*6*
                               sin(M_PI*mpSpaceMesh.ReadSpaceNode(i+1)));
}

mass.MatrixVectorMultiplier( PreviousSolution, RHS );

LHS.MatrixSolver( RHS, x );


        for (auto k: x)
            std::cout << k << ' ';

        for (auto k: x)
            myfile << k <<  " ," ;

        myfile << "\n";
        std::cout << " \n";

        for (auto k: AnalyticSolution)
            std::cout << k << ' ';
        std::cout << " \n";

        AnalyticSolution.clear();
        PreviousSolution = x;
}

        myfile.close();
}


//void HeatEquation::SetUniformSystem( double endTime, int numberOfTimeSteps, int numberOfSpaceElements,
//                            const std::string outputFileName )
//{
//    mT = endTime;
//    mm = numberOfTimeSteps;
//    mn = numberOfSpaceElements;
//    moutputFileName = outputFileName;
//}
//
//void HeatEquation::SolveUniformMesh()
//{
//
//double h = pow( mn, -1);
//
//double M_PI = 2*acos(0);
//
//double a = pow(M_PI,-2);
//
//std::vector<double> initialSpaceNodes = {0};
//std::vector<double> initialTimeNodes = {0};
//std::vector<double> initialSpaceMesh;
//std::vector<double> initialTimeMesh;
//
//
//for(int i=1; i<=mn; i++)
//    {
//        initialSpaceNodes.push_back( i*h );
//        initialSpaceMesh.push_back( initialSpaceNodes.at(i)- initialSpaceNodes.at(i-1) );
//    }
//
//for(int i=1; i<=mm; i++)
//    {
//        initialTimeNodes.push_back( i*(mT/mm) );
//        initialTimeMesh.push_back( initialTimeNodes.at(i)- initialTimeNodes.at(i-1) );
//    }
//
//
//StiffnessMatrix stiff;
//stiff.BuildStiffnessMatrix( a, initialSpaceMesh );
//
//stiff.MultiplyByScalar( initialTimeMesh.at(0) );
//
//
//
//MassMatrix mass;
//mass.BuildMassMatrix(initialSpaceMesh);
//
//TriDiagMatrix LHS;
//LHS.AddTwoMatrices( mass, stiff );
//
//
//std::vector<double> x;
//std::vector<double> AnalyticSolution;
//std::vector<double> PreviousSolution;
//std::vector<double> RHS;
//
//for (int i = 0; i<mn-1; i++)
//{
//    PreviousSolution.push_back(6*sin(M_PI*initialSpaceNodes.at(i+1)));
//}
//
//ofstream myfile;
//  myfile.open (moutputFileName);
//
//for(int j = 0; j<mm; j++)
//{
//    for (int i = 0; i<mn-1; i++)
//{
//    AnalyticSolution.push_back(exp(-initialTimeNodes.at(j+1))*6*
//                               sin(M_PI*initialSpaceNodes.at(i+1)));
//}
//
//mass.MatrixVectorMultiplier( PreviousSolution, RHS );
//
//LHS.MatrixSolver( RHS, x );
//
//
//        for (auto k: x)
//            std::cout << k << ' ';
//
//        for (auto k: x)
//            myfile << k <<  " ," ;
//
//        myfile << "\n";
//        std::cout << " \n";
//
//        for (auto k: AnalyticSolution)
//            std::cout << k << ' ';
//        std::cout << " \n";
//
//        AnalyticSolution.clear();
//        PreviousSolution = x;
//}
//
//        myfile.close();
//
//}
//
