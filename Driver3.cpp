#include <iostream>
#include <cmath>
#include "heatEquation.hpp"
#include <vector>
#include <array>
#include <math.h>

int main(int argc, char* argv[])
{

double T = 4.0;

    //n is the number of elements i.e. 1 less than the number of nodes
int n = 4;

int m = 4;

double h = pow( n, -1);

double M_PI = 2*acos(0);

double a = pow(M_PI,-2);

std::vector<double> initialSpaceNodes = {0};
std::vector<double> initialTimeNodes = {0};
std::vector<double> initialSpaceMesh;


for(int i=1; i<=n; i++)
    {
        initialSpaceNodes.push_back( i*h );
        initialTimeNodes.push_back( i*(T/m) );
        initialSpaceMesh.push_back( initialSpaceNodes.at(i)- initialSpaceNodes.at(i-1) );
    }

std::vector<double> MassDiagonal;
std::vector<double> MassLowerDiag;
std::vector<double> MassUpperDiag;

std::vector<double> StiffnessDiagonal;
std::vector<double> StiffnessLowerDiag;
std::vector<double> StiffnessUpperDiag;

HeatEquation heat;

heat.BuildStiffnessMatrix(n, a, initialSpaceMesh, StiffnessDiagonal, StiffnessLowerDiag, StiffnessUpperDiag);


        for (auto j: StiffnessDiagonal)
            std::cout << j << ' ';
        std::cout << " \n";


        for (auto j: StiffnessUpperDiag)
            std::cout << j << ' ';

        std::cout << " \n";


        for (auto j: StiffnessLowerDiag)
            std::cout << j << ' ';

std::vector<double> Diagonal;
std::vector<double> LowerDiag;
std::vector<double> UpperDiag;
std::vector<double> f;
std::vector<double> x;



f.assign( n, 1 );
Diagonal.assign( n, 3 );
UpperDiag.assign( n, 1);
UpperDiag.at( n-1 ) = 0;
LowerDiag.assign( n, 1 );
LowerDiag.at( 0 ) = 0;
x.assign( n, 0 );




}

