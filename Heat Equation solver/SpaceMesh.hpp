#ifndef SPACEMESHHEADERREF
#define SPACEMESHHEADERREF
#include <iostream>
#include <vector>
#include <array>

class SpaceMesh
{
public:

    void CopySpaceMesh (const SpaceMesh& oldSpaceMesh);
    void GenerateSpaceMesh ( std::vector<double> SpaceNodes );


    void GenerateDefaultSpaceMesh();
    void GloballyBisectSpaceMesh ();
    void InsertSpaceNode ( double xi );
    void RemoveSpaceNode (int i);
    void RefreshSpaceMesh();
    void PrintSpaceNodes ();
    void PrintSpaceMesh ();
    double ReadSpaceNode (int i);
    double ReadSpaceMesh (int i);
    int meshsize();
    int IndexAbove ( double x );


protected:

    std::vector<double> mpSpaceNodes ;
    std::vector<double> mpSpaceMesh;
    int mpmeshsize;

};


#endif // SPACEMESHHEADERREF


