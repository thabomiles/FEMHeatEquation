#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include <array>
#include "GeneralizedHeat.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <fstream>
#include <string>
#include <functional>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;


void GeneralHeat::EnergyNorm()
{
mpEnergyNorm.clear();
BuildGradientVec(mpx, mpsmesh, FEMGradient);

auto SquaredError = [this](double x)
    { return pow(GradientFunction(x) - mppde->AnalyticGradientWRTx(x, mptmesh.ReadTimeStep(mpcurrenTimeStep)), 2); };

for(int i=0; i<mpsmesh.meshsize(); i++)
{
mpEnergyNorm.push_back( gauss<double, 7>::integrate(SquaredError,
                                        mpsmesh.ReadSpaceNode(i), mpsmesh.ReadSpaceNode(i+1)) );

}
        double globalError=0;
        for(auto k: mpEnergyNorm)
            globalError = globalError + k;
        //std::cout << sqrt(globalError);
}

void GeneralHeat::AddVectors(std::vector<double>func1, std::vector<double> func2, std::vector<double>& result)
{
    if(func1.size()==func2.size())
    {
    result.clear();
     for ( int i=0; i<func1.size(); i++ )
     {
         result.push_back(func1.at(i)+func2.at(i));
     }
    }
    else
    {
        std::cout<< " your vectors are different sizes";
        std::cout<<"\n";
    }
}

void GeneralHeat::SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, APDE& apde )
{
    mpsmesh = smesh;
    mptmesh = tmesh;
    mppde = &apde;
    k_0 = mppde->k_0;
    k_L = mppde->k_L;
    g_0 = mppde->g_0;
    g_L = mppde->g_L;
    mpa = mppde->a;

    //std::cout<<apde.k_0<<"\n"<<apde.k_L<<"\n"<<mpa<<"\n"<<g_0<<"\n"<<g_L;
}

double GeneralHeat::ContinuousAnalyticSolution( double x, double t )
{
     return mppde->ContinuousAnalyticSolution( x, t );
}

void GeneralHeat::StationaryHeatEquation()
{
buildfvec( mpsmesh );

stiff.SetParameters(k_0, k_L, mpa);

stiff.BuildGeneralStiffnessMatrix ( mpsmesh );

BuiltbrVec();
AddVectors(br, f_vec, br);

stiff.MatrixSolver( br, mpx );
}

void GeneralHeat::buildfvec( SpaceMesh& a_smesh)
{
    double my_var = 0.5*a_smesh.ReadSpaceMesh(0);
    f_vec = {mppde->EllipticalRHSfunction(0)*my_var};

for(int i =1; i<mpsmesh.meshsize(); i++)
{
    my_var = 0.5*(a_smesh.ReadSpaceMesh(i)+a_smesh.ReadSpaceMesh(i-1));
    f_vec.push_back( mppde->EllipticalRHSfunction(a_smesh.ReadSpaceNode(i))*my_var );
}

my_var = 0.5*mpsmesh.ReadSpaceMesh(mpsmesh.meshsize()-1);
f_vec.push_back( mppde->EllipticalRHSfunction(a_smesh.ReadSpaceNode(a_smesh.meshsize()))*my_var );
}

void GeneralHeat::BuiltbrVec()
{
    br.assign(mpsmesh.meshsize()+1, 0);
    br.at(0) = k_0*g_0;
    br.at(mpsmesh.meshsize()) = k_L*g_L;
}

void GeneralHeat::AnalyticSolutionVec( )
{
    mpAnalyticSolution.clear();
     for (int i = 0; i<mpsmesh.meshsize()+1; i++)
{
    mpAnalyticSolution.push_back(ContinuousAnalyticSolution( mpsmesh.ReadSpaceNode(i),
                                                                        mptmesh.ReadTimeStep(mpcurrenTimeStep)));
}
}


void GeneralHeat::PrintSolution( )
{
        BuildErrorMesh();
        AnalyticSolutionVec();
//        AnalyticGradientVec();
        EnergyNorm();
        double globalError=0;
        for(auto k: mpEnergyNorm)
            globalError = globalError + k;

        std::cout << "FEM Approximation:     ";
        PrintVector(mpx);
        std::cout << "Analytic Solution:     ";
        PrintVector(mpAnalyticSolution);
        std::cout << "Error Mesh:            ";
        PrintVector(mpErrorMesh);
        std::cout << "Global Error           ";
        GlobalSpaceError();
        std::cout << "FEM Grad approx:       ";
        PrintVector(FEMGradient);
//        std::cout << "Analytic Gradient:     ";
//        PrintVector(mpTrueGradVec);
        std::cout << "Energy Error           ";
        PrintVector(mpEnergyNorm);
        std::cout << "Global Energy Error    ";
        std::cout << globalError;
        //GlobalEnergyError();
        std::cout << "\n";
}

void GeneralHeat::BuildErrorMesh()
{
mpErrorMesh.clear();

auto SquaredError = [this](double x)
    { return pow(GeneralInterpolant(x, mpx, mpsmesh ) -
    ContinuousAnalyticSolution(x, mptmesh.ReadTimeStep(mpcurrenTimeStep)), 2); };

double Q;
for(int i=0; i<mpsmesh.meshsize(); i++)
{
Q = gauss<double, 7>::integrate(SquaredError, mpsmesh.ReadSpaceNode(i), mpsmesh.ReadSpaceNode(i+1));
mpErrorMesh.push_back( Q );
}
}

double GeneralHeat::GlobalSpaceError()
{
    BuildErrorMesh();
    double globalError=0;
    for(auto k: mpErrorMesh)
        globalError = globalError + k;

    std::cout << sqrt(globalError);
    std::cout << " \n";
    return sqrt(globalError);
}

double GeneralHeat::GlobalEnergyError()
{
    EnergyNorm();
    double globalError=0;
    for(auto k: mpEnergyNorm)
        globalError = globalError + k;

    std::cout << sqrt(globalError);
    std::cout << " \n";
    return sqrt(globalError);

}

void GeneralHeat::PrintErrorMesh()
{
    BuildErrorMesh();
    PrintVector(mpErrorMesh);
}


void GeneralHeat::SolveWithBCs()
{
mpcurrenTimeStep = 0;
mpcurrentMeshIndex = 0;
AnalyticSolutionVec();
mpPreviousSolution = mpAnalyticSolution;
stiff.SetParameters(k_0, k_L, mpa);

ofstream myfile;
ofstream myfile1;
myfile.open ("solution.csv");
myfile1.open ("X.csv");



int m = mptmesh.NumberOfTimeSteps();
for(int j = 0; j<m; j++)
{
mpcurrenTimeStep = j+1;
mpcurrentMeshIndex = j;
stiff.BuildGeneralStiffnessMatrix ( mpsmesh );
stiff.MultiplyByScalar( mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );
mass.BuildGeneralMassMatrix(mpsmesh);

LHS.AddTwoMatrices( mass, stiff );

mass.MatrixVectorMultiplier( mpPreviousSolution, mpRHS );

BuiltbrVec();
VectorTimesScalar( br, mptmesh.ReadTimeMesh(mpcurrentMeshIndex) );

AddVectors( br, mpRHS, mpRHS );
LHS.MatrixSolver( mpRHS, mpx );

for (auto k: mpx)
    myfile << k << ", ";
myfile << "\n";
for (auto k: mpsmesh.mpSpaceNodes)
    myfile1 << k << ", ";
myfile1 << "\n";

mpPreviousSolution = mpx;


if (j==int(0.5*m))
{
}
}
myfile.close();
myfile1.close();

}

double GeneralHeat::GeneralInterpolant( double x, std::vector<double>& funct, SpaceMesh& relevantMesh )
{
    std::array<double, 2> firstpoint;
    std::array<double, 2> secondpoint;

    int upperindex = relevantMesh.IndexAbove( x );

    if((upperindex==1)||(upperindex==0))
    {
    firstpoint.at(0)= relevantMesh.ReadSpaceNode(0);
    firstpoint.at(1) = funct.at(0);

    secondpoint[0] = relevantMesh.ReadSpaceNode(1);
    secondpoint.at(1) = funct.at(1);
    }
    else if (upperindex == relevantMesh.meshsize())
    {
    firstpoint.at(0)= relevantMesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = funct.at(upperindex-1);

    secondpoint[0] = relevantMesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = funct.at(upperindex);
    }
    else
    {
    firstpoint.at(0)= relevantMesh.ReadSpaceNode(upperindex-1);
    firstpoint.at(1) = funct.at(upperindex-1);

    secondpoint[0] = relevantMesh.ReadSpaceNode(upperindex);
    secondpoint.at(1) = funct.at(upperindex);
    }

    long double m = (firstpoint[1]-secondpoint[1])/(firstpoint[0]-secondpoint[0]);

    return m*(x - firstpoint[0])+firstpoint[1];
}

void GeneralHeat::BuildGradientVec( std::vector<double>& funct, SpaceMesh& relevantMesh, std::vector<double>& gradvec )
{
    gradvec.clear();
    std::array<double, 2> firstpoint;
    std::array<double, 2> secondpoint;
    long double m;

    for(int i = 0; i<relevantMesh.meshsize(); i++)
    {
    firstpoint.at(0)= relevantMesh.ReadSpaceNode(i);
    firstpoint.at(1) = funct.at(i);

    secondpoint[0] = relevantMesh.ReadSpaceNode(i+1);
    secondpoint.at(1) = funct.at(i+1);

    m = (firstpoint[1]-secondpoint[1])/(firstpoint[0]-secondpoint[0]);

    gradvec.push_back(m);
    }
}

void GeneralHeat::GradientRecoveryFunction( SpaceMesh& relevantMesh,
                                             std::vector<double>& gradvec, std::vector<double>& gradrecovery )
{
    gradrecovery.clear();

    double x_0 = 0.5*(relevantMesh.ReadSpaceNode(1)+relevantMesh.ReadSpaceNode(0));
    double y_0 = gradvec.at(0);
    double x_1 = relevantMesh.ReadSpaceNode(1);
    double y_1 = 0.5*(gradvec.at(1)+gradvec.at(0));

    gradrecovery.push_back(y_0+(relevantMesh.ReadSpaceNode(0)-x_0)*(y_1-y_0)/(x_1-x_0));

    for(int i = 0; i<relevantMesh.meshsize()-1; i++)
    {
        gradrecovery.push_back(0.5*(gradvec.at(i)+gradvec.at(i+1)));
    }

    x_0 = relevantMesh.ReadSpaceNode(mpsmesh.meshsize()-1);
    y_0 = gradrecovery.back();
    x_1 = 0.5*(relevantMesh.ReadSpaceNode(relevantMesh.meshsize())+relevantMesh.ReadSpaceNode(relevantMesh.meshsize()-1));
    y_1 = gradvec.back();

    gradrecovery.push_back(y_0+(mpsmesh.ReadSpaceNode(mpsmesh.meshsize())-x_0)*(y_1-y_0)/(x_1-x_0));
}

void GeneralHeat::BuildErrorEstimate(  )
{
    ErrorEstimate.clear();

    auto GradSquaredError = [this](double x)
        { return pow(GeneralInterpolant(x, GradientRecovery, mpsmesh ) - GradientFunction(x), 2); };

    double dummy_var;
    for(int i=0; i<mpsmesh.meshsize(); i++)
    {
    ErrorEstimate.push_back( gauss<double, 7>::integrate(GradSquaredError, mpsmesh.ReadSpaceNode(i), mpsmesh.ReadSpaceNode(i+1)) );
    }
}


    //discontinuous function which throws exceptions at undefined points
double GeneralHeat::GradientFunction ( double x )
{
    int upperindex = mpsmesh.IndexAbove( x );
    if (mpsmesh.Contained(x))
    {
        std::cout<< "FEM gradient undefined at this point"<<"\n";
        return 0;
    }
    else
    {
       return FEMGradient.at(upperindex-1);
    }
}

void GeneralHeat::PrintVector( std::vector<double> aVector)
{
        for (auto k: aVector)
        std::cout << k << ", ";
        std::cout << " \n";
}

void GeneralHeat::VectorTimesScalar( std::vector<double>& func1, double scalar)
{
         for ( int i=0; i<func1.size(); i++ )
     {
         func1.at(i)= scalar*func1.at(i);
     }
}

double GeneralHeat::H_1Norm()
{
    EnergyNorm();
    BuildErrorMesh();
    double globalError=0;
    for(int i=0;i<mpErrorMesh.size(); i++)
    {
        globalError =mpEnergyNorm.at(i)+mpErrorMesh.at(i)+globalError;
    }

    std::cout << sqrt(globalError);
    std::cout << " \n";
    return sqrt(globalError);

}

void GeneralHeat::UnitTest1 ()
{
    BuildGradientVec(mpx, mpsmesh, FEMGradient);
    GradientRecoveryFunction( mpsmesh, FEMGradient, GradientRecovery );

    BuildErrorEstimate();

    double globalError=0;
    for(int i=0;i<ErrorEstimate.size(); i++)
    {
        globalError =ErrorEstimate.at(i)+globalError;
    }

    std::cout << sqrt(globalError)<<"\n";

    //PrintVector(ErrorEstimate);
}
