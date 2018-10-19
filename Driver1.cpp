#include <iostream>
#include <cmath>
#include "heatEquation.hpp"
#include <vector>
#include <array>
#include <math.h>

int main(int argc, char* argv[])
{

double T = 4.0;

int n = 3;

int m = 4;

double h = 1/n;

double M_PI = 2*acos(0);

double a = pow(M_PI,-2);

std::vector<double> initialSpaceNodes;
std::vector<double> initialTimeNodes;


std::vector<double> Diagonal = {3, 5};
std::vector<double> LowerDiag = {0, 3};
std::vector<double> UpperDiag = {4, 0};
std::vector<double> f = {2, 16};
std::vector<double> x = {0, 0};


}

