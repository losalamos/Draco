//----------------------------------*-C++-*----------------------------------//
// PoomaMesh_XYZFactory.hh
// Randy M. Roberts
// Fri Aug 20 13:33:37 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __POOMA_MT_test_PoomaMesh_XYZFactory_hh__
#define __POOMA_MT_test_PoomaMesh_XYZFactory_hh__

#include "../PoomaMesh_XYZ.hh"
#include <vector>
#include <algorithm>

namespace rtt_POOMA_MT_test
{
 
//===========================================================================//
// class PoomaMesh_XYZFactory - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class PoomaMesh_XYZFactory 
{
    // STATIC CONSTANT
    
    static const int Dimension_s = 3;
    
    // NESTED CLASSES AND TYPEDEFS

    typedef Cartesian<Dimension_s> PoomaMesh_t;
    typedef PoomaMesh_t::MeshValue_t MeshValue_t;
    
  public:

    typedef PoomaMesh_XYZ<PoomaMesh_t> MT;

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

    int numCells_m[Dimension_s];
    MeshValue_t *cellWidth_m[Dimension_s];
    e_dim_tag decomposition_m[Dimension_s];
    
  public:

    // STATIC METHODS

    static int Dimension() { return Dimension_s; }

    // CREATORS
    
    PoomaMesh_XYZFactory(const std::vector<int> &numCells_in,
			 const std::vector<std::vector<double> >
			 &cellWidth_in,
			 const std::vector<e_dim_tag> &decomposition_in
			 = std::vector<e_dim_tag>(Dimension_s, PARALLEL))
    {
	Assert(numCells_in.size() == Dimension());
	Assert(cellWidth_in.size() == Dimension());
	Assert(decomposition_in.size() == Dimension());

	std::copy(numCells_in.begin(), numCells_in.end(), numCells_m);
	for (int i=0; i<Dimension(); i++)
	{
	    Assert(cellWidth_in[i].size() == numCells_m[i]);
	    cellWidth_m[i] = new double[numCells_m[i]];
	    std::copy(cellWidth_in[i].begin(), cellWidth_in[i].end(),
		      cellWidth_m[i]);
	}
	std::copy(decomposition_in.begin(), decomposition_in.end(),
		  decomposition_m);
    }

    PoomaMesh_XYZFactory(const PoomaMesh_XYZFactory &rhs)
    {
	std::copy(rhs.numCells_m, rhs.numCells_m + Dimension(), numCells_m);
	for (int i=0; i<Dimension(); i++)
	{
	    cellWidth_m[i] = new double[numCells_m[i]];
	    std::copy(rhs.cellWidth_m[i], rhs.cellWidth_m[i] + numCells_m[i],
		      cellWidth_m[i]);
	}
	std::copy(rhs.decomposition_m, rhs.decomposition_m + Dimension(),
		  decomposition_m);
    }

    ~PoomaMesh_XYZFactory()
    {
	for (int i = 0; i < Dimension(); i++)
	{
	    delete [] cellWidth_m[i];
	}
    }

    // MANIPULATORS
    
    // ACCESSORS

    inline Product create() const;

  private:
    
    // DISSALLOWED CREATORS
    
    // DISSALLOWED MANIPULATORS
    
    PoomaMesh_XYZFactory& operator=(const PoomaMesh_XYZFactory &rhs);

    // IMPLEMENTATION
};

// Class that contains both a mesh and a field constructor.

class PoomaMesh_XYZFactory::Product
{
    rtt_dsxx::SP<MT> spmesh_m;
    MT::FieldConstructor fCtor_m;
	
  public:
	
    Product(const rtt_dsxx::SP<MT> &spmesh_in)
	: spmesh_m(spmesh_in), fCtor_m(spmesh_in)
    {
	// empty
    }
	
    MT &mesh() { return *spmesh_m; }
    MT::FieldConstructor &fieldConstructor()
    {
	return fCtor_m;
    }
};

// This inline function must be defined after nested Product class
// is defined.

inline PoomaMesh_XYZFactory::Product PoomaMesh_XYZFactory::create() const
{
    int *numCells = const_cast<int *>(numCells_m);
    MeshValue_t **cellWidth = const_cast<MeshValue_t **>(cellWidth_m);
    e_dim_tag *decomposition = const_cast<e_dim_tag *>(decomposition_m);

    return Product(rtt_dsxx::SP<MT>(new MT(numCells, cellWidth, decomposition)));
}

} // end namespace rtt_POOMA_MT_test

#endif                          // __POOMA_MT_test_PoomaMesh_XYZFactory_hh__

//---------------------------------------------------------------------------//
//                              end of POOMA_MT/test/PoomaMesh_XYZFactory.hh
//---------------------------------------------------------------------------//
