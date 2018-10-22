#include <iostream>
#include <vector>
#include <array>
#include "TriDiagMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <cmath>

void StiffnessMatrix::BuildStiffnessMatrix( double a, std::vector<double> SpaceMesh )
{
mpn = SpaceMesh.size()-1;
mpDiagonal = { a*pow(SpaceMesh.at(0), -1) + a*pow(SpaceMesh.at(1), -1) };
mpLowerDiag = {0};
mpUpperDiag = { -1*a*pow( SpaceMesh.at(1), -1 ) };

for(int i=1; i<mpn-1; i++)
    {
        mpDiagonal.push_back( a*pow(SpaceMesh.at(i), -1) + a*pow(SpaceMesh.at(i+1), -1) );
        mpLowerDiag.push_back( -1*a*pow(SpaceMesh.at(i),-1)  );
        mpUpperDiag.push_back( -1*a*pow(SpaceMesh.at(i+1), -1) );
    }

mpDiagonal.push_back( a*pow(SpaceMesh.at(mpn-1), -1) + a*pow(SpaceMesh.at(mpn-1), -1) );
mpLowerDiag.push_back( -1*a*pow(SpaceMesh.at(mpn-1),-1) );
mpUpperDiag.push_back( 0 );
}
