#ifndef LU_CONFIG_H
#define LU_CONFIG_H

#include "pup.h"

#ifndef STOP_AFTER
    #define STOP_AFTER -1
#endif

class LUConfig {
    public:
        LUConfig():
            blockSize(0), numBlocks(0), pivotBatchSize(0),
            mappingScheme(0), peTileRows(0), peTileCols(0), peTileRotate(0),
            memThreshold(0)
        {}


        void pup(PUP::er &p)
        {
            p | blockSize;
            p | numBlocks;
            p | pivotBatchSize;

            p | mappingScheme;
            p | peTileRows;
            p | peTileCols;
            p | peTileRotate;

            p | memThreshold;
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
        int mappingScheme;
        /// The num of PEs in the PEtile
        int peTileRows;
        /// The num of rows in the PEtile
        int peTileCols;
        /// The amount of rotation in the PEtile
        int peTileRotate;
        /// The max memory (in MB) that the adaptive RTS should limit the app to
        int memThreshold;
};

#endif // LU_CONFIG_H

