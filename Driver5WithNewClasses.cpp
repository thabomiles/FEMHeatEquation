#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "heatEquation.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"

int main(int argc, char* argv[])
{

double T = 1.0;

    //n is the number of elements i.e. 1 less than the number of nodes
int n = 4;

    //number of time steps i.e. 1 less than the number of nodes
int m = 4;

double h = pow( n, -1);

double M_PI = 2*acos(0);

double a = pow(M_PI,-2);

std::vector<double> initialSpaceNodes = {0};
std::vector<double> initialTimeNodes = {0};
std::vector<double> initialSpaceMesh;
std::vector<double> initialTimeMesh;


for(int i=1; i<=n; i++)
    {
        initialSpaceNodes.push_back( i*h );
        initialSpaceMesh.push_back( initialSpaceNodes.at(i)- initialSpaceNodes.at(i-1) );
    }

for(int i=1; i<=m; i++)
    {
        initialTimeNodes.push_back( i*(T/m) );
        initialTimeMesh.push_back( initialTimeNodes.at(i)- initialTimeNodes.at(i-1) );
    }


StiffnessMatrix stiff;
stiff.BuildStiffnessMatrix( a, initialSpaceMesh );

stiff.MultiplyByScalar( initialTimeMesh.at(0) );



MassMatrix mass;
mass.BuildMassMatrix(initialSpaceMesh);

TriDiagMatrix LHS;
LHS.AddTwoMatrices( mass, stiff );


std::vector<double> x;
std::vector<double> AnalyticSolution;
std::vector<double> PreviousSolution;
std::vector<double> RHS;

for (int i = 0; i<n-1; i++)
{
    PreviousSolution.push_back(6*sin(M_PI*initialSpaceNodes.at(i+1)));
}

for(int j = 0; j<m; j++)
{
    for (int i = 0; i<n-1; i++)
{
    AnalyticSolution.push_back(exp(-initialTimeNodes.at(j+1))*6*
                               sin(M_PI*initialSpaceNodes.at(i+1)));
}

mass.MatrixVectorMultiplier( PreviousSolution, RHS );

LHS.MatrixSolver( RHS, x );


        for (auto j: x)
            std::cout << j << ' ';

        std::cout << " \n";


        for (auto j: AnalyticSolution)
            std::cout << j << ' ';
        std::cout << " \n";

        AnalyticSolution.clear();
        PreviousSolution = x;
}



}

