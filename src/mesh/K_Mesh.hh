//----------------------------------*-C++-*----------------------------------//
// K_Mesh.hh
// Geoffrey Furnish
// Wed Jun 11 00:14:08 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __mesh_K_Mesh_hh__
#define __mesh_K_Mesh_hh__

#include "ds++/SP.hh"
#include "ds++/Mat.hh"

#include "K_Mesh_DB.hh"

#include "Pooma.h"

//===========================================================================//
// class K_Mesh - A 1-d structured mesh

// 
//===========================================================================//

template<class Geometry>
class K_Mesh : private Geometry,
	       private K_Mesh_DB
{
  public:
    Index C, F;
    FieldLayout<1> clayout, flayout;

  public:

// cell centered scalar field

    class ccsf : public Field<double, 1>
    {
	SP< K_Mesh<Geometry> > m;

	ccsf( K_Mesh<Geometry> *pm );
      public:
	ccsf( SP< K_Mesh<Geometry> >& m_ );

	friend class K_Mesh<Geometry>;
    };

// face centered scalar field

    class fcsf : public Field<double, 1>
    {
	SP< K_Mesh<Geometry> > m;

	fcsf( K_Mesh<Geometry> *pm );
      public:
	fcsf( SP< K_Mesh<Geometry> >& m_ );
	fcsf( FieldLayout<1>& fl ) : Field<double,1>( fl ) {}

	friend class K_Mesh<Geometry>;
    };

  private:
    ccsf xc, Vc;		// cell center coord, cell volume.
    fcsf xf, Af;		// face coord, face area.

  public:
    K_Mesh( const K_Mesh_DB& kdb );
};

#endif                          // __mesh_K_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/K_Mesh.hh
//---------------------------------------------------------------------------//
