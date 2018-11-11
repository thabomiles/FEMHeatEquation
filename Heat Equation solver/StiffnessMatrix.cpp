#include <iostream>
#include <vector>
#include <array>
#include "TriDiagMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <cmath>

void StiffnessMatrix::BuildStiffnessMatrix( SpaceMesh smesh )
{
mpDiagonal.clear();
mpLowerDiag.clear();
mpUpperDiag.clear();

mpn = smesh.meshsize()-1;
mpDiagonal = { a*pow(smesh.ReadSpaceMesh(0), -1) + a*pow(smesh.ReadSpaceMesh(1), -1) };
mpLowerDiag = {0};
mpUpperDiag = { -1*a*pow( smesh.ReadSpaceMesh(1), -1 ) };

for(int i=1; i<mpn-1; i++)
    {
        mpDiagonal.push_back( a*pow(smesh.ReadSpaceMesh(i), -1) + a*pow(smesh.ReadSpaceMesh(i+1), -1) );
        mpLowerDiag.push_back( -1*a*pow(smesh.ReadSpaceMesh(i),-1)  );
        mpUpperDiag.push_back( -1*a*pow(smesh.ReadSpaceMesh(i+1), -1) );
    }

mpDiagonal.push_back( a*pow(smesh.ReadSpaceMesh(mpn-1), -1) + a*pow(smesh.ReadSpaceMesh(mpn), -1) );
mpLowerDiag.push_back( -1*a*pow(smesh.ReadSpaceMesh(mpn-1),-1) );
mpUpperDiag.push_back( 0 );
}

void StiffnessMatrix::BuildGeneralStiffnessMatrix ( SpaceMesh smesh )
{
mpDiagonal.clear();
mpLowerDiag.clear();
mpUpperDiag.clear();

mpn = smesh.meshsize()+1;
mpDiagonal = { a*pow(smesh.ReadSpaceMesh(0), -1)+mpk_0 };
mpLowerDiag = {0};
mpUpperDiag = { -1*a*pow( smesh.ReadSpaceMesh(0), -1 ) };

for(int i=0; i<mpn-2; i++)
    {
        mpDiagonal.push_back( a*pow(smesh.ReadSpaceMesh(i), -1) + a*pow(smesh.ReadSpaceMesh(i+1), -1) );
        mpLowerDiag.push_back( -1*a*pow(smesh.ReadSpaceMesh(i),-1)  );
        mpUpperDiag.push_back( -1*a*pow(smesh.ReadSpaceMesh(i+1), -1) );
    }

mpDiagonal.push_back( a*pow(smesh.ReadSpaceMesh(mpn-2), -1)+mpk_L );
mpLowerDiag.push_back( -1*a*pow(smesh.ReadSpaceMesh(mpn-2),-1) );
mpUpperDiag.push_back( 0 );

}

void StiffnessMatrix::SetParameters (double k_0, double k_L)
{
    mpk_0=k_0;
    mpk_L=k_L;
}

