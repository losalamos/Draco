//----------------------------------*-C++-*----------------------------------//
// MeshXYZConnFacesAroundVertices.hh
// Randy M. Roberts
// Wed Apr 28 10:19:19 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __mesh_MeshXYZConnFacesAroundVertices_hh__
#define __mesh_MeshXYZConnFacesAroundVertices_hh__

#include "./MeshXYZConnFAVIterators.hh"

namespace rtt_mesh
{
 
//===========================================================================//
// class MeshXYZConnFacesAroundVertices - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class FaceField>
class MeshXYZConnFacesAroundVertices 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename FaceField::value_type FaceValType;
    typedef NSVertex::Vertex<FaceField> VertexType;

  public:
    
    typedef VertexType value_type;

    typedef MeshXYZConnFAVIterator<FaceField> iterator;
    typedef const_MeshXYZConnFAVIterator<FaceField> const_iterator;
    
  private:

    // DATA

    FaceField &field;
    int numCells;

  public:

    // CREATORS
    
    MeshXYZConnFacesAroundVertices(FaceField &field_)
	: field(field_)
    {
	numCells = field.size() / NSVertex::NUMFACES_PER_CELL;
    }
    MeshXYZConnFacesAroundVertices(const MeshXYZConnFacesAroundVertices &rhs)
	: field(rhs.field)
    {
	numCells = field.size() / NSVertex::NUMFACES_PER_CELL;
    }
    ~MeshXYZConnFacesAroundVertices() { /* empty */ }

    // MANIPULATORS
    
    MeshXYZConnFacesAroundVertices&
    operator=(const MeshXYZConnFacesAroundVertices &rhs);

    // ACCESSORS

    iterator begin()
    {
	return iterator(field, 0, 0);
    }
    iterator end()
    {
	return iterator(field, numCells, 0);
    }

    const_iterator begin() const
    {
	return iterator(field, 0, 0);
    }
    const_iterator end() const
    {
	return iterator(field, numCells, 0);
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_mesh

#endif                          // __mesh_MeshXYZConnFacesAroundVertices_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/MeshXYZConnFacesAroundVertices.hh
//---------------------------------------------------------------------------//
