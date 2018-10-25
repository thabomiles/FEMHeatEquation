#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include "SpaceMesh.hpp"

void SpaceMesh::GenerateSpaceMesh( std::vector<double> SpaceNodes )
{
    mpSpaceNodes = SpaceNodes;
    RefreshSpaceMesh();
    mpmeshsize = mpSpaceMesh.size();
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
    std::cout << k << ' ';
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
    return mpSpaceMesh.at(i);
}

int SpaceMesh::meshsize()
{
    return mpSpaceMesh.size();
}


