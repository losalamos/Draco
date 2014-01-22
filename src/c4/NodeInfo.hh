//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   NodeInfo.hh
 * \author Geoffrey Furnish
 * \date   Tue Jan 17 10:13:47 1995
 * \brief  Class to hold parallel configuration information.
 * \note   Copyright (C) 1995-2014 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_NodeInfo_hh__
#define __c4_NodeInfo_hh__

#include "C4_Functions.hh"

namespace rtt_c4
{

//===========================================================================//
/*! \class NodeInfo - Parallel configuration information
 *
 * This class contains information about the configuration of a parallel
 * multicomputer.  User objects may inherit from this in order to learn where
 * they fit into the total scheme of things.
 */
class DLL_PUBLIC NodeInfo
{

  public:
    int node;
    int nodes;
    int lastnode;

    NodeInfo(void)
        : node( rtt_c4::node() ),
          nodes( rtt_c4::nodes() ),
          lastnode( nodes-1 )
    { /* empty */  };
    virtual ~NodeInfo(void) {/*empty*/};
};

} // end of rtt_c4

#endif // __c4_NodeInfo_hh__

//---------------------------------------------------------------------------//
// end of c4/NodeInfo.hh
//---------------------------------------------------------------------------//
