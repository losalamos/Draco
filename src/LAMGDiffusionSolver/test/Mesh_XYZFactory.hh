//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZFactory.hh
// Randy M. Roberts
// Fri Aug 20 13:33:37 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_test_Mesh_XYZFactory_hh__
#define __LAMGDiffusionSolver_test_Mesh_XYZFactory_hh__

#include "mesh/Mesh_XYZ.hh"
#include <string>

namespace rtt_LAMGDiffusionSolver_test
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

    // Nested Product Class to be defined below.
    
    class Product;

    // These classes are tags required by some clients.

    struct Structured { /* empty */ };
    struct UnStructured { /* empty */ };

    // This Structuring typedef is used for the client to switch,
    // at compile time, between which methods are instantiated and invoked.
    //
    // This particular mesh factory produces a "Structured" mesh.
    
    typedef Structured Structuring;

  private:
    
    // DATA

    const Mesh_DB mdb_m;
    
  public:

    // STATIC METHODS

    static int Dimension() { return 3; }
    static std::string name() { return "Mesh_XYZFactory"; }

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

    inline Product create() const;

  private:
    
    // DISSALLOWED CREATORS
    
    // Mesh_XYZFactory(const Mesh_XYZFactory &rhs);

    // DISSALLOWED MANIPULATORS
    
    Mesh_XYZFactory& operator=(const Mesh_XYZFactory &rhs);

    // IMPLEMENTATION
};

// Class that contains both a mesh and a field constructor.

class Mesh_XYZFactory::Product
{
    rtt_dsxx::SP<Mesh_XYZ> spmesh_m;
    Mesh_XYZ::FieldConstructor fCtor_m;
	
  public:
	
    Product(const rtt_dsxx::SP<Mesh_XYZ> &spmesh_in)
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

// This inline function must be defined after nested Product class
// is defined.

inline Mesh_XYZFactory::Product Mesh_XYZFactory::create() const
{
    return Product(rtt_dsxx::SP<Mesh_XYZ>(new Mesh_XYZ(mdb_m)));
}

} // end namespace rtt_LAMGDiffusionSolver_test

#endif // __LAMGDiffusionSolver_test_Mesh_XYZFactory_hh__

//---------------------------------------------------------------------------//
// end of LAMGDiffusionSolver/test/Mesh_XYZFactory.hh
//---------------------------------------------------------------------------//
