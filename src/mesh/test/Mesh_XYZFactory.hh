//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZFactory.hh
// Randy M. Roberts
// Fri Aug 20 13:33:37 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __mesh_test_Mesh_XYZFactory_hh__
#define __mesh_test_Mesh_XYZFactory_hh__

#include "../Mesh_XYZ.hh"

namespace rtt_mesh_test
{
 
//===========================================================================//
// class Mesh_XYZFactory - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Mesh_XYZFactory 
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef Mesh_XYZ MT;
    
    class Product
    {
	dsxx::SP<Mesh_XYZ> spmesh_m;
	Mesh_XYZ::FieldConstructor fCtor_m;
	
      public:
	
	Product(const dsxx::SP<Mesh_XYZ> &spmesh_in)
	    : spmesh_m(spmesh_in), fCtor_m(spmesh_in)
	{
	    // empty
	}
	
	Mesh_XYZ &mesh() { return *spmesh_m; }
	Mesh_XYZ::FieldConstructor &fieldConstructor()
	{
	    return fCtor_m;
	}
    };

  private:
    
    // DATA

    const Mesh_DB &mdb_m;
    
  public:

    // CREATORS
    
    Mesh_XYZFactory(const Mesh_DB &mdb_in)
	: mdb_m(mdb_in)
    {
	// empty
    }

    ~Mesh_XYZFactory()
    {
	// empty
    }

    // MANIPULATORS
    
    // ACCESSORS

    Product create() const
    {
	return Product(dsxx::SP<Mesh_XYZ>(new Mesh_XYZ(mdb_m)));
    }

  private:
    
    // DISSALLOWED CREATORS
    
    Mesh_XYZFactory(const Mesh_XYZFactory &rhs);

    // DISSALLOWED MANIPULATORS
    
    Mesh_XYZFactory& operator=(const Mesh_XYZFactory &rhs);

    // IMPLEMENTATION
};

} // end namespace rtt_mesh_test

#endif                          // __mesh_test_Mesh_XYZFactory_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/test/Mesh_XYZFactory.hh
//---------------------------------------------------------------------------//
