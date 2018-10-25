#include "GeneralMesh.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>

void GeneralMesh::GeneralSpaceMesh( std::vector<double> SpaceNodes )
{
    mpSpaceNodes = SpaceNodes;
    RefreshSpaceMesh();
}

void GeneralMesh::GeneralTimeMesh ( std::vector<double> TimeNodes )
{
    mpTimeNodes = TimeNodes;
    RefreshTimeMesh();
}

void GeneralMesh::GloballyBisectSpaceMesh ()
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
void GeneralMesh::GloballyBisectTimeMesh ()
{
    double gridPoint;
    for (int i=0; i<mpTimeMesh.size(); i++)
    {
        gridPoint = mpTimeNodes.at(i)+0.5*mpTimeMesh.at(i);
        mpTimeNodes.push_back(gridPoint);
    }
    sort(mpTimeNodes.begin(), mpTimeNodes.end());

    RefreshTimeMesh();
}
void GeneralMesh::InsertSpaceNode ( double xi )
{
    mpSpaceNodes.push_back( xi );
    sort(mpSpaceNodes.begin(), mpSpaceNodes.end());

    RefreshSpaceMesh();
}

void GeneralMesh::InsertTimeNode ( double ti )
{
    mpTimeNodes.push_back( ti );
    sort(mpTimeNodes.begin(), mpTimeNodes.end());

    RefreshTimeMesh();
}

void GeneralMesh::RefreshSpaceMesh()
{
    mpSpaceMesh.clear();
    for (int i = 0; i<mpSpaceNodes.size()-1; i++)
    mpSpaceMesh.push_back( mpSpaceNodes.at(i+1)- mpSpaceNodes.at(i) );
}

void GeneralMesh::RefreshTimeMesh()
{
     mpTimeMesh.clear();
    for (int i = 0; i<mpTimeNodes.size()-1; i++)
    mpTimeMesh.push_back( mpTimeNodes.at(i+1)- mpTimeNodes.at(i) );
}


void GeneralMesh::PrintSpaceNodes ()
{
    for (auto k: mpSpaceNodes)
    std::cout << k << ' ';
    std::cout << " \n";
}
void GeneralMesh::PrintTimeNodes ()
{
    for (auto k: mpTimeNodes)
    std::cout << k << ' ';
    std::cout << " \n";
}

void GeneralMesh::PrintSpaceMesh ()
{
    for (auto k: mpSpaceMesh)
    std::cout << k << ' ';
    std::cout << " \n";
}

void GeneralMesh::PrintTimeMesh ()
{
    for (auto k: mpTimeMesh)
    std::cout << k << ' ';
    std::cout << " \n";
}


