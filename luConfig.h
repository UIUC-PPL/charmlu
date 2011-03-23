#ifndef LU_CONFIG_H
#define LU_CONFIG_H

#include "pup.h"
#include <charm++.h>

#ifndef STOP_AFTER
    #define STOP_AFTER -1
#endif

class LUConfig {
    public:
        LUConfig():
            blockSize(0), numBlocks(0), pivotBatchSize(0),
            peTileRows(0), peTileCols(0), peTileRotate(0), peTileStride(0),
            memThreshold(0),
            tracePeriodFraction(0.25), numStepsToTrace(3),
            numTimesToTrace(1 + 1.0/tracePeriodFraction)
        {}


        void pup(PUP::er &p)
        {
            p | blockSize;
            p | numBlocks;
            p | pivotBatchSize;

            p | map;
            p | peTileRows;
            p | peTileCols;
            p | peTileRotate;
            p | peTileStride;

            p | memThreshold;

            p | tracePeriodFraction;
            p | numStepsToTrace;
            p | numTimesToTrace;
        }


        /// The size of the matrix
        int matrixSize;
        /// The size of each square block of the matrix
        int blockSize;
        /// The number of blocks comprising the whole square matrix
        int numBlocks;
        /// The number of pivots to agglomerate into one batch
        int pivotBatchSize;
        /// The mapping logic that has to be used
        CkGroupID map;
        /// The num of PEs in the PEtile
        int peTileRows;
        /// The num of rows in the PEtile
        int peTileCols;
        /// The amount of rotation in the PEtile
        int peTileRotate;
        /// The stride between PEs on a column
        int peTileStride;
        /// The max memory (in MB) that the adaptive RTS should limit the app to
        int memThreshold;

        /// The period (as a fraction of numBlocks) with which to toggle tracing
        double tracePeriodFraction;
        /// The number of active panels to trace each time
        int numStepsToTrace;
        /// The number of times when tracing will be active (computed from above two)
        int numTimesToTrace;
};

#endif // LU_CONFIG_H

