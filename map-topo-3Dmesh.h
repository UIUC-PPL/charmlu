#include "peMeshDims.h"
#include "TopoManager.h"

#ifndef MAP_TOPO_3D_MESH_H
#define MAP_TOPO_3D_MESH_H

/// A global topomgr object for everyone to use
extern TopoManager luTopoMgr;

/**
 * A mapping for the LU chare array tuned for machines with 3D meshes / torii
 *
 * Mimics a 2D tile mapping by creating a 2D tile out of a 3D mesh of
 * processors.  Each column of the chare array (panel) is mapped onto a
 * complete YZ plane of nodes in the 3D mesh The number of cores / node used
 * for each panel can also be configured.
 * Succesive panels are be mapped onto
 *    - succesive YZ processor slices along the X direction of the mesh
 *    - wrap around onto a different set of cores on the node after all the slices are used once
 *
 * Thus, some cores in each YZ slice of nodes host a panel. This becomes the
 * number of rows in the 2D peTile.  Each slice needs to host a number  of
 * panels to use all the cores.  The number of YZ slices (X dimension of the 3D
 * mesh) x the number of panels per slice is the number of columns in the 2D
 * peTile.
 *
 * Motivations:
 * - Minimize interference  of other msgs with the active panel
 *      Active panel msgs (small) are likely contending with block msgs
 *      (large). Avoiding this may help accelerate the active panel.  This can
 *      be ahieved by mapping a panel along a plane of the mesh.  Communication
 *      within the panel should be restricted to this plane.  However, other
 *      panels pulling or pushing data should mostly use links in perpendicular
 *      directions.
 *
 *      The above assumption will be perfectly valid if we use all the PEs on a
 *      node to host a panel. If this is not the case, then there can be non
 *      active panel PEs on a node that interact with other similar nodes on
 *      the same plane of the machine. However, hopefully this can be
 *      controlled and minimized.
 *
 * - Allow balancing memory and network demands of the active panel
 *      Some experiments seem to indicate that it is possible to use too many PEs
 *      to host a panel, which stresses the memory resources on that node when the
 *      panel becomes active. The extent of the effect is different on different machines.
 *      Hence, its desirable to use only some PEs on a node to host a panel.
 *
 *      At the same time, we also do not want to put succesive panels on the
 *      same nodes, as this can still cause interference between the pulls on
 *      the L blocks from the previous active panel with the next active panel.
 *      There is also the minor interference of the original active panel with
 *      the agglomerated pivoting on the next panel (which is also quite
 *      important). Hence its preferable to put succesive panels on different
 *      sets of nodes.
 *
 *  This set of requirements has a possible solution in the naive, slice-based
 *  mapping below.
 */
class LUMapTopo: public LUMap
{
    public:
        ///
        LUMapTopo(const int _numBlocks, PEMeshDims panelPEmesh):
            numBlocks(_numBlocks),
            allPEdims(luTopoMgr.getDimNX(), luTopoMgr.getDimNY(), luTopoMgr.getDimNZ(), luTopoMgr.getDimNT()),
            activePanelPEdims(panelPEmesh)

        {
            if (_numBlocks <= 0)
                CkAbort("numBlocks < 0!!");

            if ( !isFeasible() )
                CkAbort("Active panels cannot be mapped onto the specified sub-mesh. Check your inputs to the mapping scheme.");

            if (CkMyPe() == 0)
            {
                CkPrintf("\tMachine Topology: %d x %d x %d x %d\n", allPEdims.x, allPEdims.y, allPEdims.z, allPEdims.t);
                CkPrintf("\tPanel PE Mesh: %d x %d x %d x %d\n", activePanelPEdims.x, activePanelPEdims.y, activePanelPEdims.z, activePanelPEdims.t);
            }
            // Compute the number of panels (chare columns) per plane of the 3d pe mesh
            numPanelsPerMeshPlane = allPEdims.t / activePanelPEdims.t;
            // Compute the total number of panels (chare columns) that will fit onto the
            // whole 3d pe mesh. This constitutes the number of columns in the pe tile
            numPanelsInTile = allPEdims.x * numPanelsPerMeshPlane;
        }



        /// Return the PE where this index should reside
        inline int procNum(int arrayHdl, const CkArrayIndex &arrIdx)
        {
            const int *idx = arrIdx.data();
            int x, y, z, t;
            // Assume that an active panel is mapped onto a complete YZ plane
            int linearizedYZidx = idx[0]/activePanelPEdims.t;
            int ta = idx[0]  % activePanelPEdims.t;
            int tb = (idx[1]%numPanelsInTile)  / allPEdims.x;

            x  = idx[1]  % allPEdims.x;
            y  = linearizedYZidx / activePanelPEdims.z;
            z  = linearizedYZidx % activePanelPEdims.z;
            t  = ta * numPanelsPerMeshPlane + tb;

            CkAssert( (x >= 0 && x < allPEdims.x) &&
                      (y >= 0 && y < allPEdims.y) &&
                      (z >= 0 && z < allPEdims.z) &&
                      (t >= 0 && t < allPEdims.t)
                    );

            return luTopoMgr.coordinatesToRank(x, y, z, t);
        }



        /// Check if its feasible to map the active panel onto the specified sub-mesh
        inline bool isFeasible() const
        {
            // Ensure this is a valid mesh size
            if (!activePanelPEdims.isValid())
                return false;

            // Ensure that user has asked to map the active panel onto
            // at most a 2D mesh of nodes. Not more (3D)
            if ( (activePanelPEdims.x > 1) &&
                 (activePanelPEdims.y > 1) &&
                 (activePanelPEdims.z > 1)
               )
                return false;

            // Ensure that each node can host an integer number of active panels
            // This massively simplifies load balance and mapping
            if (allPEdims.t % activePanelPEdims.t != 0)
                return false;

            // To start with,
            // Assume that the active panel is always mapped onto the XY plane
            // Hence, check if one YZ plane can host an integral number of active panels
            // in each dimension
            if ( (allPEdims.y % activePanelPEdims.y != 0) ||
                 (allPEdims.z % activePanelPEdims.z != 0)
               )
                return false;

            // Assume that each panel can be load balanced onto a complete YZ plane
            if (numBlocks % (allPEdims.y * allPEdims.z) == 0 && CkMyPe() == 0)
                CkPrintf("WARNING: mapping will not yield load balance. Hence uneven mem usage\n");
            // Assume that the panels can be evenly distributed across the whole PE space
            if (numBlocks % ( allPEdims.x * (allPEdims.t / activePanelPEdims.t) ) == 0 && CkMyPe() == 0)
                CkPrintf("WARNING: mapping will not yield load balance. Hence uneven mem usage\n");
            return true;
        }



        ///
        inline int pesInPanel(CkIndex2D index)
        {
            return activePanelPEdims.x * activePanelPEdims.y * activePanelPEdims.z * activePanelPEdims.t; 
        }

        ///
        inline std::string desc() { return "LU array mapping for 3D meshes and torii"; }

    private:
        /// Machine dimensions
        PEMeshDims allPEdims;
        /// Dimensions of partition onto which an active panel is mapped
        PEMeshDims activePanelPEdims;
        // The number of panels (chare columns) per plane of the 3d pe mesh
        int numPanelsPerMeshPlane;
        /// The number of panels that will fit onto the whole 3d pe mesh
        int numPanelsInTile;
        /// The total number of blocks in the A matrix
        int numBlocks;
};

// determine X, Y, Z and T dims
// compute XY, XZ, XT, YZ, YT, ZT
// Find factors of T. Ta and Tb such that Ta . Tb == T
// Find longest dimension (Da) of X,Y,Z
// PE tile dimensions should be (Dl.Ta, Db
// We dont want active panel spread along T dimension
// peRows has to be one of these numbers or a factor or multiple for performance
#endif // MAP_TOPO_3D_MESH_H

