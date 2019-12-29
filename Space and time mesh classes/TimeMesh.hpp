#ifndef TIMEMESHHEADERREF
#define TIMEMESHHEADERREF
#include <iostream>
#include <vector>
#include <array>

class TimeMesh
{
public:

    void CopyTimeMesh (const TimeMesh& oldTimeMesh);
    void GenerateTimeMesh ( std::vector<double> TimeNodes );

        //takes mTimeStep which is one less than the number of
        //nodes and a Final time and produces the nodes i.e. m+1 equal steps
    void GenerateUniformTimeMesh ( int mTimeSteps, double TFinalTime );
    void BisectInterval (int lowerIndex, int upperIndex);
    void GloballyBisectTimeMesh ();
    void InsertTimeNode ( double ti );
    void RefreshTimeMesh();
    void PrintTimeNodes ();
    void PrintTimeMesh ();
    double ReadTimeStep (int i);
    double ReadTimeMesh (int i);
    int NumberOfTimeSteps();
    std::vector<double> mpTimeNodes;



protected:

    int mpnumberOfTimeSteps;

    std::vector<double> mpTimeMesh;

};


#endif // TIMEMESHHEADERREF


