#include "TopoManager.h"
#include "pup.h"

#ifndef PE_MESH_DIMS_H
#define PE_MESH_DIMS_H

/// A global topomgr object for everyone to use
extern TopoManager luTopoMgr;

/// Utility struct that represents the dimensions of a 3D processor mesh as found in machines like BGL, BGP and Cray XT5 and earlier
class PEMeshDims
{
    public:
        int x, y, z, t;

        PEMeshDims():
            x(0), y(0), z(0), t(0)
        {}

        PEMeshDims(const int _x, const int _y, const int _z, const int _t):
            x(_x), y(_y), z(_z), t(_t)
        {}

        void pup(PUP::er &p)
        {
            p|x;
            p|y;
            p|z;
            p|t;
        }

        inline bool isValid() const
        {
            if (x <=0 || y <= 0 || z <= 0 || t <= 0)
                return false;

            const int dX = luTopoMgr.getDimNX();
            const int dY = luTopoMgr.getDimNY();
            const int dZ = luTopoMgr.getDimNZ();
            const int dT = luTopoMgr.getDimNT();

            if (x > dX || y > dY || z > dZ || t > dT)
                return false;

            return true;
        }
};

#endif // PE_MESH_DIMS_H

