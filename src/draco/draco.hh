//----------------------------------*-C++-*----------------------------------//
// draco.hh
// Geoffrey Furnish
// Thu Jul 31 15:55:22 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __draco_draco_hh__
#define __draco_draco_hh__

#include "c4/NodeInfo.hh"

#include "geom/Cartesian1d.hh"
#include "mesh/K_Mesh.hh"

#include "draco/draconml.hh"

#include "ds++/SP.hh"

//===========================================================================//
// class draco - 

// 
//===========================================================================//

class draco : public NodeInfo, public draconml
{
    SP< K_Mesh<Cartesian1d> > mesh;
  public:
    draco();

    SP< K_Mesh<Cartesian1d> > get_mesh() { return mesh; }
};

#endif                          // __draco_draco_hh__

//---------------------------------------------------------------------------//
//                              end of draco/draco.hh
//---------------------------------------------------------------------------//
