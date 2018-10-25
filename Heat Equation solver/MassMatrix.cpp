#include <iostream>
#include <vector>
#include <array>
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include <cmath>
#include "SpaceMesh.hpp"

void MassMatrix::BuildMassMatrix ( SpaceMesh smesh )
{
mpn = smesh.meshsize()-1;
mpDiagonal = {(smesh.ReadSpaceMesh(0)+smesh.ReadSpaceMesh(1))*pow(3,-1)};
mpLowerDiag = {0};
mpUpperDiag = {(smesh.ReadSpaceMesh(1)*pow(6,-1))};
//mpDiagonal = {(SpaceMesh.at(0)+SpaceMesh.at(1))*pow(3,-1)};
//mpLowerDiag = {0};
//mpUpperDiag = {(SpaceMesh.at(1)*pow(6,-1))};

for(int i=1; i<mpn-1; i++)
    {
        mpDiagonal.push_back( (smesh.ReadSpaceMesh(i)+smesh.ReadSpaceMesh(i+1))*pow(3,-1) );
        mpLowerDiag.push_back( smesh.ReadSpaceMesh(i)*pow(6,-1) );
        mpUpperDiag.push_back( (smesh.ReadSpaceMesh(i+1)*pow(6,-1)) );
    }

mpDiagonal.push_back( (smesh.ReadSpaceMesh(mpn-1)+smesh.ReadSpaceMesh(mpn))*pow(3,-1) );
mpLowerDiag.push_back( smesh.ReadSpaceMesh(mpn-1)*pow(6,-1) );
mpUpperDiag.push_back( 0 );
         }

