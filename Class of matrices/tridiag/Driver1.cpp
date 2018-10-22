#include <iostream>
#include <cmath>
#include "TriDiagMatrix.hpp"
#include <vector>
#include <array>
#include <math.h>

int main(int argc, char* argv[])
{

int n = 4;



double M_PI = 2*acos(0);

double a = pow(M_PI,-2);

std::vector<double> UpperDiag = {4, 0, 0};
std::vector<double> Diagonal = {3, 5, 3};
std::vector<double> LowerDiag = {0, 3, 0};
std::vector<double> f = {2, 16, 1};
std::vector<double> x = {0, 0};

TriDiagMatrix matrix;

matrix.SetMatrix(n-1, Diagonal, LowerDiag, UpperDiag);

matrix.PrintMatrix();

matrix.MatrixVectorMultiplier(f,x);

        for (auto j: x)
        std::cout << j << ' ';
        std::cout << " \n";

}

