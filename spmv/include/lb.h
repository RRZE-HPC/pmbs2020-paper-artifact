#ifndef RACE_LB_H
#define RACE_LB_H

#include "type.h"
#include "error.h"

enum neighbour_t{
    acquire,
    give
};

class Stat{
    private:
        int *arr;
        int len;
        std::vector<int> scale;

    public:
        Stat(int *arr_, int len_, int numPartitions_, std::vector<int> scale_);
        ~Stat();
        int numPartitions;
        void calculate();
        double *mean;
        double *var;
        double *acquireWeight;
        double *giveWeight;
        double *weight;
};

class LB{
    private:
        int* levelPtr;
        int* subLevelPtr;

        /////////// Overall layout ///////////////////////////
        //  |   zonePtr[0]  subZonePtr[0]                   //
        //  |               subZonePtr[1]                   //
        //  |               ...                             //
        //  |               subZonePtr[numBlocks[0]-1]      //
        //  |   zonePtr[1]  subZonePtr[numBlocks[0]]        //
        //  |               subZonePtr[numBlocks[0]+1]      //
        //  |               ...                             //
        //  |               subZonePtr[numBlocks[1]-1]      //
        //  |   ...         ...                             //
        //////////////////////////////////////////////////////
        int* subZonePtr;//required for blocking
        int* numBlocks;
        int* zonePtr;

        double* effRatio;
        std::vector<int> scale;

        double* targetData; 

        //gap to be left between 2 levels
        int minGap;

        int maxThreads;
        int nThreads;
        int currLvlThreads;
        double efficiency;
        int blockPerThread;
        int totalBlocks;
        int totalSubBlocks;

        RACE::LBTarget lbTarget;

        void calcChunkSum_general(int *arr, int *levelPtr_, int len, bool forceRow=false);
        void calcChunkSum(int *arr, bool forceRow=false);
        void calcSubChunkSum(int *arr, bool forceRow=false);
        void calcEffectiveRatio();
        void calcZonePtr(int base);
        void calcSubZonePtr(int base);
        int  findNeighbour(const Stat &stats, neighbour_t type);
        void moveOneStep(int toIdx, int fromIdx);
        void splitZones();

    public:
        LB(int nThreads_, double efficiency_, LevelData* levelData_, RACE::dist dist_, RACE::d2Method d2Type=RACE::TWO_BLOCK, RACE::LBTarget lbTarget=RACE::NNZ); //constructor
        ~LB(); //destructor

        int getMaxThreads();
        int getNumThreads();
        void getZonePtr(int **zonePtr_, int *len, int base=0);
        void getSubZonePtr(int **subZonePtr_, int *len, int base=0);
        void getNumBlocks(int **numBlocks_, int *len);
        void getNxtLvlThreads(int **nxtLvlThreads, int *len);
        int getCurrLvlNThreads();
        RACE_error balance();
};


#endif
