#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <string>
#include "TimeMesh.hpp"

void TimeMesh::CopyTimeMesh (const TimeMesh& oldTimeMesh)
{
    mpTimeNodes = oldTimeMesh.mpTimeNodes;
    RefreshTimeMesh();
}

void TimeMesh::BisectInterval (int lowerIndex, int upperIndex)
{
    if ((upperIndex-lowerIndex==1)&&(upperIndex<mpTimeNodes.size()))
    {
    double midpoint = 0.5*(mpTimeNodes.at(upperIndex)+mpTimeNodes.at(lowerIndex));
    mpTimeNodes.push_back( midpoint );
    sort(mpTimeNodes.begin(), mpTimeNodes.end());
    }
    else
    {
    std::cout<< "this is not an interval";
    std::cout<< "\n";
    }

    RefreshTimeMesh();
}

void TimeMesh::GenerateTimeMesh ( std::vector<double> TimeNodes )
{
    mpTimeNodes = TimeNodes;
    RefreshTimeMesh();
}

void TimeMesh::GenerateUniformTimeMesh ( int mTimeSteps, double TFinalTime )
{
    mpTimeNodes = {0};
    for(int i=1; i<=mTimeSteps; i++)
    {
        mpTimeNodes.push_back( i*(TFinalTime/mTimeSteps) );
    }
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

int TimeMesh::NumberOfTimeSteps()
{
    return mpTimeMesh.size();
}

