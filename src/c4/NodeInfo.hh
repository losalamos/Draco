//----------------------------------*-C++-*----------------------------------//
// NodeInfo.hh
// Geoffrey Furnish
// Tue Jan 17 10:13:47 1995
//---------------------------------------------------------------------------//
// @> Class to hold parallel configuration information.
//---------------------------------------------------------------------------//

#ifndef __c4_NodeInfo_hh__
#define __c4_NodeInfo_hh__

#include "c4/config.hh"

//===========================================================================//
// class NodeInfo - Parallel configuration information

// This class contains information about the configuration of a parallel
// multicomputer.  User objects may inherit from this in order to learn where
// they fit into the total scheme of things.
//===========================================================================//

class NodeInfo {

  public:
    int node;
    int nodes;
    int group;

    int lastnode;

    NodeInfo();
};

#endif                          // __c4_NodeInfo_hh__

//---------------------------------------------------------------------------//
//                              end of c4/NodeInfo.hh
//---------------------------------------------------------------------------//
