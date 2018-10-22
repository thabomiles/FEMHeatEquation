#include <iostream>
#include <vector>
#include <array>
#include "TriDiagMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <cmath>

void StiffnessMatrix::BuildStiffnessMatrix( double a, std::vector<double> SpaceMesh )
{
mpn = SpaceMesh.size();
mpDiagonal = { a*pow(SpaceMesh.at(0), -1) + a*pow(SpaceMesh.at(1), -1) };
mpLowerDiag = {0};
mpUpperDiag = { -1*a*pow( SpaceMesh.at(1), -1 ) };

for(int i=1; i<mpn-2; i++)
    {
        mpDiagonal.push_back( a*pow(SpaceMesh.at(i), -1) + a*pow(SpaceMesh.at(i+1), -1) );
        mpLowerDiag.push_back( -1*a*pow(SpaceMesh.at(i),-1)  );
        mpUpperDiag.push_back( -1*a*pow(SpaceMesh.at(i+1), -1) );
    }

mpDiagonal.push_back( a*pow(SpaceMesh.at(mpn-2), -1) + a*pow(SpaceMesh.at(mpn-3), -1) );
mpLowerDiag.push_back( -1*a*pow(SpaceMesh.at(mpn-2),-1) );
mpUpperDiag.push_back( 0 );
}
