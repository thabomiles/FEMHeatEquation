#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include "TimeMesh.hpp"

void TimeMesh::GenerateTimeMesh ( std::vector<double> TimeNodes )
{
    mpTimeNodes = TimeNodes;
    RefreshTimeMesh();
}

void TimeMesh::GloballyBisectTimeMesh ()
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

void TimeMesh::InsertTimeNode ( double ti )
{
    mpTimeNodes.push_back( ti );
    sort(mpTimeNodes.begin(), mpTimeNodes.end());

    RefreshTimeMesh();
}

void TimeMesh::RefreshTimeMesh()
{
     mpTimeMesh.clear();
    for (int i = 0; i<mpTimeNodes.size()-1; i++)
    mpTimeMesh.push_back( mpTimeNodes.at(i+1)- mpTimeNodes.at(i) );
}


void TimeMesh::PrintTimeNodes ()
{
    for (auto k: mpTimeNodes)
    std::cout << k << ' ';
    std::cout << " \n";
}

void TimeMesh::PrintTimeMesh ()
{
    for (auto k: mpTimeMesh)
    std::cout << k << ' ';
    std::cout << " \n";
}

double TimeMesh::ReadTimeStep(int i)
{
    return mpTimeNodes.at(i);
}


double TimeMesh::ReadTimeMesh (int i)
{
    return mpTimeMesh.at(i);
}



