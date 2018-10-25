#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "heatEquation.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include "GeneralMesh.hpp"
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char* argv[])
{

double T = 1.0;

    //n is the number of elements i.e. 1 less than the number of nodes

    //number of time steps i.e. 1 less than the number of nodes
int m = 16;

//HeatEquation heat = HeatEquation( T, m, n, "soultion1.csv");
//
//heat.SolveUniformMesh();


//double h = pow( n, -1);

double M_PI = 2*acos(0);

double a = pow(M_PI,-2);

std::vector<double> initialSpaceNodes = {0, 0.25, 0.5, 1};
int n = initialSpaceNodes.size()-1;
std::vector<double> initialTimeNodes = {0};
std::vector<double> initialSpaceMesh;
std::vector<double> initialTimeMesh;


//for(int i=1; i<=n; i++)
//    {
//        initialSpaceNodes.push_back( i*h );
//        initialSpaceMesh.push_back( initialSpaceNodes.at(i)- initialSpaceNodes.at(i-1) );
//    }

for(int i=1; i<=m; i++)
    {
        initialTimeNodes.push_back( i*(T/m) );
        initialTimeMesh.push_back( initialTimeNodes.at(i)- initialTimeNodes.at(i-1) );
    }

SpaceMesh smesh;
smesh.GenerateSpaceMesh(initialSpaceNodes);
smesh.GloballyBisectSpaceMesh();
smesh.GloballyBisectSpaceMesh();

n=smesh.meshsize();



TimeMesh tmesh;
tmesh.GenerateTimeMesh( initialTimeNodes );
smesh.PrintSpaceNodes();
smesh.PrintSpaceMesh();
std::cout << n << ' ';

StiffnessMatrix stiff;
stiff.BuildStiffnessMatrix( a, smesh );
//stiff.PrintMatrix();


stiff.MultiplyByScalar( tmesh.ReadTimeMesh(0) );
//stiff.PrintMatrix();





MassMatrix mass;
mass.BuildMassMatrix(smesh);
//mass.PrintMatrix();
TriDiagMatrix LHS;
LHS.AddTwoMatrices( mass, stiff );


std::vector<double> x;
std::vector<double> AnalyticSolution;
std::vector<double> PreviousSolution;
std::vector<double> RHS;

for (int i = 0; i<n-1; i++)
{
    PreviousSolution.push_back(6*sin(M_PI*smesh.ReadSpaceNode(i+1)));
}

ofstream myfile;
  myfile.open ("soultion1.txt");

for(int j = 0; j<m; j++)
{
    for (int i = 0; i<n-1; i++)
{
    AnalyticSolution.push_back(exp(-tmesh.ReadTimeStep(j+1))*6*
                               sin(M_PI*smesh.ReadSpaceNode(i+1)));
}

mass.MatrixVectorMultiplier( PreviousSolution, RHS );

LHS.MatrixSolver( RHS, x );


        for (auto k: x)
            std::cout << k << ' ';

        for (auto k: x)
            myfile << k <<  " ," ;

        myfile << "\n";
        std::cout << " \n";




        for (auto p: AnalyticSolution)
            std::cout << p << ' ';
        std::cout << " \n";




        AnalyticSolution.clear();
        PreviousSolution = x;
}

        myfile.close();

}

