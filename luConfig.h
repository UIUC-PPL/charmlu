#include "pup.h"

#ifndef LU_CONFIG_H
#define LU_CONFIG_H

class LUConfig {
    public:
        LUConfig():
            blockSize(0), numBlocks(0), pivotBatchSize(0),
            memThreshold(0)
        {}


        void pup(PUP::er &p)
        {
            p | blockSize;
            p | numBlocks;
            p | pivotBatchSize;

            p | memThreshold;
        }



        /// The size of each square block of the matrix
        int blockSize;
        /// The number of blocks comprising the whole square matrix
        int numBlocks;
        /// The number of pivots to agglomerate into one batch
        int pivotBatchSize;
        /// The max memory (in MB) that the adaptive RTS should limit the app to
        int memThreshold;
};

#endif // LU_CONFIG_H

