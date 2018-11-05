#include <iostream>
#include <cmath>
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include <vector>
#include <array>
#include <math.h>
#include "StiffnessMatrix.hpp"

void printx (std::vector<double> x)
{
    for (auto i: x)
        std::cout<<i<< ", ";
    std::cout<< "\n";
}

int main(int argc, char* argv[])
{

int n = 4;
double h = pow(n, -1);

double M_PI = 2*acos(0);

double a = pow(M_PI,-2);

std::vector<double> UpperDiag = {7, -10, 0};
std::vector<double> Diagonal = {3, 7, +10};
std::vector<double> LowerDiag = {0, -3, 1};
std::vector<double> f = {15, 13, 3};
std::vector<double> x = {1, 2, 3};
std::vector<double> SpaceMesh;
SpaceMesh.assign(n, h);








MassMatrix mass1;

mass1.BuildMassMatrix( SpaceMesh );

mass1.MultiplyByScalar( 24 );


mass1.PrintMatrix();

mass1.MatrixSolver(f, x);
printx(f);
mass1.PrintMatrix();


printx(x);








}



