#include <iostream>
#include <vector>
#include <array>
#include "GeneralMesh.hpp"

int main(int argc, char* argv[])
{

double T = 4;

int n = 3;

int m = 4;

double h = 1/n;

std::vector<double> SpaceNodes = {0, 0.5, 0.75, 1};
std::vector<double> TimeNodes = {0, 0.1, 0.2, 0.5, 1};

GeneralMesh  mesh;
mesh.GeneralSpaceMesh( SpaceNodes );
mesh.GeneralTimeMesh( TimeNodes );

mesh.PrintSpaceNodes();
mesh.PrintSpaceMesh();

mesh.GloballyBisectSpaceMesh();

//mesh.PrintSpaceNodes();
//mesh.PrintSpaceMesh();
//mesh.InsertSpaceNode( 0.66 );
//mesh.PrintSpaceNodes();

mesh.PrintTimeNodes();
mesh.PrintTimeMesh();
mesh.GloballyBisectTimeMesh();
mesh.PrintTimeNodes();
mesh.PrintTimeMesh();
mesh.InsertTimeNode( 0.3333 );
mesh.PrintTimeNodes();
mesh.PrintTimeMesh();
}
