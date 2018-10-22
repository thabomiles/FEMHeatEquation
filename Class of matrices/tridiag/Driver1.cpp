#include <iostream>
#include <cmath>
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include <vector>
#include <array>
#include <math.h>
#include "StiffnessMatrix.hpp"

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
std::vector<double> x = {0, 0, 0};
std::vector<double> SpaceMesh;
SpaceMesh.assign(n, h);



//TriDiagMatrix matrix;
//
//matrix.SetMatrix(n-1, Diagonal, LowerDiag, UpperDiag);
//
//matrix.PrintMatrix();
//
//matrix.MatrixVectorMultiplier(f,x);
//
//        for (auto j: x)
//        std::cout << j << ' ';
//        std::cout << " \n";

//MassMatrix mass;
//
//mass.BuildMassMatrix( SpaceMesh );
//
//mass.PrintMatrix();

//StiffnessMatrix stiffness;
//
//stiffness.BuildStiffnessMatrix( a, SpaceMesh);
//
//stiffness.PrintMatrix();

//TriDiagMatrix matrix;
//
//matrix.SetMatrix( Diagonal, LowerDiag, UpperDiag);
//
//matrix.PrintMatrix();
//
//matrix.MultiplyByScalar(2);
//
//matrix.PrintMatrix();

MassMatrix mass1, mass2;
TriDiagMatrix matrix;

mass1.BuildMassMatrix( SpaceMesh );
mass2.BuildMassMatrix( SpaceMesh);


mass1.MultiplyByScalar( 24 );
mass2.MultiplyByScalar( 24 );

mass1.PrintMatrix();

matrix.AddTwoMatrices(mass1, mass2);

matrix.PrintMatrix();







}

