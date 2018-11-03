#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include "SpaceMesh.hpp"

void SpaceMesh::CopySpaceMesh (const SpaceMesh& oldSpaceMesh)
{

    mpSpaceNodes = oldSpaceMesh.mpSpaceNodes;
    RefreshSpaceMesh();
    mpmeshsize = mpSpaceMesh.size();
}

void SpaceMesh::BisectInterval (int lowerIndex, int upperIndex)
{
    if ((upperIndex-lowerIndex==1)&&(upperIndex<mpSpaceNodes.size()))
    {
    double midpoint = 0.5*(mpSpaceNodes.at(upperIndex)+mpSpaceNodes.at(lowerIndex));
    mpSpaceNodes.push_back( midpoint );
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());
    }
    else
    {
    std::cout<< "this is not an interval";
    std::cout<< "\n";
    }

    RefreshSpaceMesh();
}

void SpaceMesh::GenerateSpaceMesh( std::vector<double> SpaceNodes )
{
    mpSpaceNodes = SpaceNodes;
    RefreshSpaceMesh();
    mpmeshsize = mpSpaceMesh.size();
}

void SpaceMesh::GenerateDefaultSpaceMesh()
{
    mpSpaceNodes = {0, 0.25, 0.5, 1};
    RefreshSpaceMesh();
}

void SpaceMesh::GloballyBisectSpaceMesh ()
{
    double gridPoint;
    for (int i=0; i<mpSpaceMesh.size(); i++)
    {
        gridPoint = mpSpaceNodes.at(i)+0.5*mpSpaceMesh.at(i);
        mpSpaceNodes.push_back(gridPoint);
    }
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());

    RefreshSpaceMesh();
}

void SpaceMesh::InsertSpaceNode ( double xi )
{
    mpSpaceNodes.push_back( xi );
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());

    RefreshSpaceMesh();
}

void SpaceMesh::RefreshSpaceMesh()
{
    mpSpaceMesh.clear();
    for (int i = 0; i<mpSpaceNodes.size()-1; i++)
    mpSpaceMesh.push_back( mpSpaceNodes.at(i+1)- mpSpaceNodes.at(i) );
}

void SpaceMesh::PrintSpaceNodes ()
{
    for (auto k: mpSpaceNodes)
    std::cout << k << ", ";
    std::cout << " \n";
}

void SpaceMesh::PrintSpaceMesh ()
{
    for (auto k: mpSpaceMesh)
    std::cout << k << ' ';
    std::cout << " \n";
}

double SpaceMesh::ReadSpaceNode (int i)
{
    return mpSpaceNodes.at(i);
}


double SpaceMesh::ReadSpaceMesh (int i)
{
    return mpSpaceMesh[i];
}

int SpaceMesh::meshsize()
{
    return mpSpaceMesh.size();
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
    RefreshSpaceMesh();
}



//double h = pow( n, -1);
//for(int i=1; i<=n; i++)
//    {
//        initialSpaceNodes.push_back( i*h );
//        initialSpaceMesh.push_back( initialSpaceNodes.at(i)- initialSpaceNodes.at(i-1) );
//    }


