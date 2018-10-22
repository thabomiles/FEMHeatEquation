#include <iostream>
#include <vector>
#include <array>
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include <cmath>

void MassMatrix::BuildMassMatrix ( std::vector<double> SpaceMesh )
{
mpn = SpaceMesh.size();
mpDiagonal = {(SpaceMesh.at(0)+SpaceMesh.at(1))*pow(3,-1)};
mpLowerDiag = {0};
mpUpperDiag = {(SpaceMesh.at(1)*pow(6,-1))};

for(int i=1; i<mpn-2; i++)
    {
        mpDiagonal.push_back( (SpaceMesh.at(i)+SpaceMesh.at(i+1))*pow(3,-1) );
        mpLowerDiag.push_back( SpaceMesh.at(i)*pow(6,-1) );
        mpUpperDiag.push_back( (SpaceMesh.at(i+1)*pow(6,-1)) );
    }

mpDiagonal.push_back( (SpaceMesh.at(mpn-2)+SpaceMesh.at(mpn-3))*pow(3,-1) );
mpLowerDiag.push_back( SpaceMesh.at(mpn-2)*pow(6,-1) );
mpUpperDiag.push_back( 0 );
         }

