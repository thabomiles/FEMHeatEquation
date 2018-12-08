#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include "SpaceMesh.hpp"
#include <cmath>
#include <set>

void SpaceMesh::Range( double lowerlimit, double upperlimit, std::vector<double>& Nodes)
{
    Nodes.clear();
    for(auto i:mpSpaceNodes)
        if(lowerlimit<=i&&i<=upperlimit)
        {
            Nodes.push_back(i);
        }
}

void SpaceMesh::ReadLastSpaceNode()
{
    mpSpaceNodes.back();
}


void SpaceMesh::CommonMesh( SpaceMesh& firstmesh, SpaceMesh& secondmesh )
{
    mpSpaceNodes=firstmesh.mpSpaceNodes;
    mpSpaceNodes.insert( mpSpaceNodes.end(), secondmesh.mpSpaceNodes.begin(), secondmesh.mpSpaceNodes.end() );
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());
    mpSpaceNodes.erase( unique( mpSpaceNodes.begin(), mpSpaceNodes.end() ), mpSpaceNodes.end() );
}

double SpaceMesh::GeneralTestFunctions(int nodeIndex, double x)
{
    if(nodeIndex==0)
    {
        //std::cout<< " half hat1 "<< "\n";
        return -pow(ReadSpaceMesh(nodeIndex),-1)*(x-mpSpaceNodes.at(nodeIndex))+1;
    }
    else if(nodeIndex==meshsize())
    {
        //std::cout<< " half hat2 "<< "\n";
        return pow(ReadSpaceMesh(nodeIndex-1),-1)*(x-mpSpaceNodes.at(nodeIndex))+1;
    }
    else if (x<=mpSpaceNodes[nodeIndex])
    {
        //std::cout<< " Lower "<< "\n";
        return pow(ReadSpaceMesh(nodeIndex-1),-1)*(x-mpSpaceNodes.at(nodeIndex))+1;
    }
    else if (x>mpSpaceNodes[nodeIndex])
    {
        //std::cout<< " Upper "<< "\n";
        return -pow(ReadSpaceMesh(nodeIndex),-1)*(x-mpSpaceNodes.at(nodeIndex))+1;
    }
}

double SpaceMesh::TestFunctions( int nodeIndex, double x)
{
        if((nodeIndex==0)||(nodeIndex==meshsize()))
    {
        std::cout<< " This function has 0 dirichlet boundaries "<< "\n";
    }
    else if (x<mpSpaceNodes[nodeIndex]-ReadSpaceMesh(nodeIndex-1)||x>mpSpaceNodes[nodeIndex]+ReadSpaceMesh(nodeIndex))
    {
        std::cout<< " You are out of your interval "<< "\n";
        return 0;
    }
    else if (x<=mpSpaceNodes[nodeIndex])
    {
        //std::cout<< " Lower "<< "\n";
        return pow(ReadSpaceMesh(nodeIndex-1),-1)*(x-mpSpaceNodes.at(nodeIndex))+1;
    }
    else if (x>mpSpaceNodes[nodeIndex])
    {
        //std::cout<< " Upper "<< "\n";
        return -pow(ReadSpaceMesh(nodeIndex),-1)*(x-mpSpaceNodes.at(nodeIndex))+1;
    }
}

void SpaceMesh::GenerateUniformMesh(double boundary, int numberofnodes)
{
    for(int i=0; i<numberofnodes; i++)
    {
        mpSpaceNodes.push_back( i*boundary*pow( numberofnodes-1, -1) );
    }
}

void SpaceMesh::CopySpaceMesh (const SpaceMesh& oldSpaceMesh)
{

    mpSpaceNodes = oldSpaceMesh.mpSpaceNodes;
    mpmeshsize = mpSpaceNodes.size()-1;
}

void SpaceMesh::InsertArray ( std::vector<double> RefinementNodes )
{
    for( auto i: RefinementNodes)
        mpSpaceNodes.push_back(i);
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());
    mpSpaceNodes.erase( unique( mpSpaceNodes.begin(), mpSpaceNodes.end() ), mpSpaceNodes.end() );
}

void SpaceMesh::BisectIntervals (std::vector<int> &intervalsForBisection)
{
    for(auto i: intervalsForBisection)
        mpSpaceNodes.push_back(0.5*(mpSpaceNodes.at(i)+mpSpaceNodes.at(i+1)));
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());
}

void SpaceMesh::CoarsenIntervals (std::vector<int> &intervalsForCoarsening)
{
    for (int i = intervalsForCoarsening.size()-1; i>=0; i-- )
        std::remove(mpSpaceNodes.begin(), mpSpaceNodes.end(),mpSpaceNodes.at(intervalsForCoarsening.at(i)));

    int vecSize = mpSpaceNodes.size()-intervalsForCoarsening.size();
    mpSpaceNodes.resize(vecSize);
}

void SpaceMesh::GenerateSpaceMesh( std::vector<double> SpaceNodes )
{
    mpSpaceNodes = SpaceNodes;
    mpmeshsize = meshsize();
}

void SpaceMesh::GenerateDefaultSpaceMesh()
{
    mpSpaceNodes = {0.0, 0.25, 0.5, 1.0};
}

void SpaceMesh::GloballyBisectSpaceMesh ()
{
    double mesh_size = mpSpaceNodes.size()-1;
    double gridPoint;
    for (int i=0; i<mesh_size; i++)
    {
        gridPoint = mpSpaceNodes.at(i)+0.5*ReadSpaceMesh(i);
        mpSpaceNodes.push_back(gridPoint);
    }
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());
}

void SpaceMesh::InsertSpaceNode ( double xi )
{
    mpSpaceNodes.push_back( xi );
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());
}

//void SpaceMesh::RefreshSpaceMesh()
//{
//    mpSpaceMesh.clear();
//    for (int i = 0; i<mpSpaceNodes.size()-1; i++)
//    mpSpaceMesh.push_back( mpSpaceNodes.at(i+1)- mpSpaceNodes.at(i) );
//    mpmeshsize = mpSpaceMesh.size();
//}

void SpaceMesh::PrintSpaceNodes ()
{
    for (auto k: mpSpaceNodes)
    std::cout << k << ", ";
    std::cout << " \n";
}

void SpaceMesh::PrintSpaceMesh ()
{
    for (int i = 0; i<mpSpaceNodes.size()-1; i++)
    {
        std::cout << mpSpaceNodes.at(i+1)-mpSpaceNodes.at(i) << ", ";
    }
    std::cout << " \n";
}

double SpaceMesh::ReadSpaceNode (int i)
{
    return mpSpaceNodes.at(i);
}


double SpaceMesh::ReadSpaceMesh (int i)
{
    return mpSpaceNodes.at(i+1)-mpSpaceNodes.at(i);
}

int SpaceMesh::meshsize()
{
    return mpSpaceNodes.size()-1;
}

int SpaceMesh::IndexAbove ( double x )
{
    auto it = lower_bound(mpSpaceNodes.begin(),mpSpaceNodes.end(), x);

   return distance(mpSpaceNodes.begin(), it);

}

void SpaceMesh::RemoveSpaceNode (int i)
{
    if (i==0||i==meshsize())
    {
        std::cout <<"\n";
        std::cout<< meshsize();
        std::cout <<"\n";
        std::cout<< "you cannot remove the first of last node";
        std::cout <<"\n";
    }
    else
    {
    mpSpaceNodes.erase(mpSpaceNodes.begin() + i);
    }
}

bool SpaceMesh::Contained (double my_var)
{
    bool found = (std::find(mpSpaceNodes.begin(), mpSpaceNodes.end(), my_var) != mpSpaceNodes.end());
}


