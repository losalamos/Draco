//----------------------------------*-C++-*----------------------------------//
// TraceInfo.hh
// Geoffrey Furnish
// 16 June 1994
//---------------------------------------------------------------------------//
// @> A class for holding info needed during Tcl trace processing.
//---------------------------------------------------------------------------//

#ifndef __nml_TraceInfo_hh__
#define __nml_TraceInfo_hh__

#include "Item.hh"
#include "Block.hh"

//===========================================================================//
// class NML_TraceInfo

// This class is just a holding point for keeping track of all the
// different little snipets of information which are useful for Tcl
// trace processing.
//===========================================================================//

class NML_TraceInfo {
    NML_Block *pb;
    NML_Item *pi;

  public:
    NML_TraceInfo( NML_Block *_pb, NML_Item *_pi );
    NML_Block *Get_pb() { return pb; }
    NML_Item  *Get_pi() { return pi; }
};

#endif                          // __nml_TraceInfo_hh__

//---------------------------------------------------------------------------//
//                              end of TraceInfo.hh
//---------------------------------------------------------------------------//
