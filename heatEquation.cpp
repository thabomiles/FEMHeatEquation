#include <iostream>
#include <cmath>
#include "heatEquation.hpp"
#include <vector>
#include <array>
#include <math.h>

void HeatEquation::Solve()
{

}

void HeatEquation::TridiagonalMatrixSolver( int n,
         std::vector<double> Diagonal, std::vector<double> LowerDiag,
         std::vector<double> UpperDiag, std::vector<double> f,
         std::vector<double> &x )
    {
                  //scaling matrix entries
    for(int i=1; i<n; i++)
    {
        Diagonal.at(i) = Diagonal.at(i)-(UpperDiag.at(i-1)*LowerDiag.at(i)/Diagonal.at(i-1));
        f.at(i) = f.at(i)-(f.at(i-1)*LowerDiag.at(i)/Diagonal.at(i-1));
    }

        //solving through elimination
     x.at(n-1)=f.at(n-1)/Diagonal.at(n-1);
     for(int i=n-2; i>=0; i=i-1)
     {
         x.at(i) = (f.at(i) - UpperDiag.at(i)*x.at(i+1))/(Diagonal.at(i));
     }
    }


void HeatEquation::SetSystem( std::vector<double> Diagonal,
            std::vector<double> LowerDiag, std::vector<double> UpperDiag,
            std::vector<double> f, std::vector<double> x )
    {
        mCurrentDiagonal = Diagonal;
        mCurrentLowerDiag = LowerDiag;
        mCurrentUpperDiag = UpperDiag;
        mCurrentf = f;
        mCurrentSolution = x;

    }


void HeatEquation::SetVariables( int n, int m, double T, double a )
    {
        mn = n;
        mm = m;
        mT = T;
        ma = a;

    }

void HeatEquation::SetSpaceTimeMesh( std::vector<double> currentSpaceNodes,
                           std::vector<double> currentTimeNodes )
    {
        mCurrentSpaceNodes = currentSpaceNodes;
        mCurrentTimeNodes = currentTimeNodes;

    }

void HeatEquation::MatrixVectorMultiplier ( int n, std::vector<double> Diagonal,
        std::vector<double> LowerDiag, std::vector<double> UpperDiag,
        std::vector<double> f, std::vector<double> &ProductVector )
    {
        ProductVector.assign( n, 0 );

        ProductVector.at(0) = Diagonal.at(0)*f.at(0)+UpperDiag.at(0)*f.at(1);
        ProductVector.at(n-1) = LowerDiag.at(n-1)*f.at(n-2)+Diagonal.at(n-1)*f.at(n-1);

            for (int j = 1; j<n-1; j++)
            {
                ProductVector.at(j) = LowerDiag.at(j)*f.at(j)+Diagonal.at(j)*f.at(j)+UpperDiag.at(j)*f.at(j+1);
            }

    }
    //Builds a Mass matrix of dimensions (n-1)x(n-1)
void HeatEquation::BuildMassMatrix ( int n, std::vector<double> SpaceMesh,
         std::vector<double> &MassDiagonal, std::vector<double> &MassLowerDiag,
         std::vector<double> &MassUpperDiag  )
         {

MassDiagonal = {(SpaceMesh.at(0)+SpaceMesh.at(1))*pow(3,-1)};
MassLowerDiag = {0};
MassUpperDiag = {(SpaceMesh.at(1)*pow(6,-1))};

for(int i=1; i<n-2; i++)
    {
        MassDiagonal.push_back( (SpaceMesh.at(i)+SpaceMesh.at(i+1))*pow(3,-1) );
        MassLowerDiag.push_back( SpaceMesh.at(i)*pow(6,-1) );
        MassUpperDiag.push_back( (SpaceMesh.at(i+1)*pow(6,-1)) );
    }

MassDiagonal.push_back( (SpaceMesh.at(n-2)+SpaceMesh.at(n-3))*pow(3,-1) );
MassLowerDiag.push_back( SpaceMesh.at(n-2)*pow(6,-1) );
MassUpperDiag.push_back( 0 );
         }

       //builds a stiffness matrix of size (n-1)*(n-1) for 0 boundary conditions
void HeatEquation::BuildStiffnessMatrix ( int n, double a, std::vector<double> SpaceMesh,
         std::vector<double> &StiffnessDiagonal, std::vector<double> &StiffnessLowerDiag,
         std::vector<double> &StiffnessUpperDiag  )
         {

StiffnessDiagonal = { a*pow(SpaceMesh.at(0), -1) + a*pow(SpaceMesh.at(1), -1) };
StiffnessLowerDiag = {0};
StiffnessUpperDiag = { -1*a*pow( SpaceMesh.at(1), -1 ) };

for(int i=1; i<n-2; i++)
    {
        StiffnessDiagonal.push_back( a*pow(SpaceMesh.at(i), -1) + a*pow(SpaceMesh.at(i+1), -1) );
        StiffnessLowerDiag.push_back( -1*a*pow(SpaceMesh.at(i),-1)  );
        StiffnessUpperDiag.push_back( -1*a*pow(SpaceMesh.at(i+1), -1) );
    }

StiffnessDiagonal.push_back( a*pow(SpaceMesh.at(n-2), -1) + a*pow(SpaceMesh.at(n-3), -1) );
StiffnessLowerDiag.push_back( -1*a*pow(SpaceMesh.at(n-2),-1) );
StiffnessUpperDiag.push_back( 0 );
         }

//void HeatEquation::TimeStepper
