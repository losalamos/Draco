//----------------------------------*-C++-*----------------------------------//
// MeshXYZConnFAVIterators.hh
// Randy M. Roberts
// Thu Apr 29 09:44:54 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __mesh_MeshXYZConnFAVIterators_hh__
#define __mesh_MeshXYZConnFAVIterators_hh__

#include <iterator>
#include "ds++/Assert.hh"

namespace rtt_mesh
{

namespace NSVertex
{

template<class T>
class Vertex;

enum { NUMFACES_PER_VERTX = 3, NUMVERTICES_PER_CELL = 8,
       NUMFACES_PER_CELL = 6 };
}   

template<class FaceField>
class MeshXYZConnFAVIterator_base
{
    typedef typename FaceField::value_type FaceValType;
    typedef NSVertex::Vertex<FaceField> VertexType;

    FaceField &field;
    int cellIndex;
    int vertexIndex;

    mutable VertexType vertex;

  protected:

    MeshXYZConnFAVIterator_base(FaceField &field_, int cellIndex_,
				int vertexIndex_)
	: field(field_), cellIndex(cellIndex_), vertexIndex(vertexIndex_),
	  vertex(field, cellIndex, vertexIndex)
    {
	// empty
    }

    MeshXYZConnFAVIterator_base(const MeshXYZConnFAVIterator_base &vit)
	: field(vit.field), cellIndex(vit.cellIndex),
	  vertexIndex(vit.vertexIndex),
	  vertex(field, cellIndex, vertexIndex)
    {
	// empty
    }

    virtual ~MeshXYZConnFAVIterator_base()
    {
	// empty
    }

  public:
    
    bool operator==(const MeshXYZConnFAVIterator_base &rhs) const
    {
	return vertexIndex == rhs.vertexIndex && cellIndex == rhs.cellIndex;
    }
    bool operator!=(const MeshXYZConnFAVIterator_base &rhs) const
    {
	return !operator==(rhs);
    }

  protected:
    
    template<class Reference>
    Reference operator*() const
    {
	Assert (vertexIndex >= 0);
	Assert (vertexIndex < NSVertex::NUMVERTICES_PER_CELL);
	Assert (cellIndex >= 0);
	// number of cells = field.size() / NUMFACES_PER_CELL.
	Assert (cellIndex < field.size() / NSVertex::NUMFACES_PER_CELL);
	return vertex;
    }

    template<class Pointer>
    Pointer operator->() const
    {
	Assert (vertexIndex >= 0);
	Assert (vertexIndex < NSVertex::NUMVERTICES_PER_CELL);
	Assert (cellIndex >= 0);
	// number of cells = field.size() / NUMFACES_PER_CELL.
	Assert (cellIndex < field.size() / NUMFACES_PER_CELL);
	return &vertex;
    }

    void incr()
    {
	vertexIndex++;
	if (vertexIndex == NSVertex::NUMVERTICES_PER_CELL)
	{
	    cellIndex++;
	    vertexIndex = 0;
	}
	vertex.assign(cellIndex, vertexIndex);
    }
};

// Forward Declaration

template<class FaceField>
class const_MeshXYZConnFAVIterator;

template<class FaceField>
class MeshXYZConnFAVIterator
    : public std::iterator<std::forward_iterator_tag, int,
      NSVertex::Vertex<FaceField>, NSVertex::Vertex<FaceField> *,
      NSVertex::Vertex<FaceField> &>,
      public MeshXYZConnFAVIterator_base<FaceField>
{
    friend class const_MeshXYZConnFAVIterator<FaceField>;

#if 0
    typedef std::iterator<std::forward_iterator_tag, int,
      NSVertex::Vertex<FaceField>, NSVertex::Vertex<FaceField> *,
	NSVertex::Vertex<FaceField> &> Iterator;

    typedef typename Iterator::reference reference;
    typedef typename Iterator::pointer pointer;
#endif
    
  public:

    MeshXYZConnFAVIterator(FaceField &field_, int cellIndex_, int vertexIndex_)
	: MeshXYZConnFAVIterator_base<FaceField>(field_, cellIndex_,
						 vertexIndex_)
    {
	// empty
    }

    MeshXYZConnFAVIterator(const MeshXYZConnFAVIterator &it_)
	: MeshXYZConnFAVIterator_base<FaceField>(it_) { /* empty */ }
	
    MeshXYZConnFAVIterator &operator++()
    {
	incr(); return *this;
    }
    MeshXYZConnFAVIterator operator++(int)
    {
	MeshXYZConnFAVIterator tmp = *this; incr(); return tmp;
    }

    typename MeshXYZConnFAVIterator::reference operator*()
    { return MeshXYZConnFAVIterator_base<FaceField>::operator*<reference>(); }
    typename MeshXYZConnFAVIterator::pointer operator->()
    { return MeshXYZConnFAVIterator_base<FaceField>::operator-><pointer>(); }
};
    
template<class FaceField>
class const_MeshXYZConnFAVIterator
    : public std::iterator<std::forward_iterator_tag, int,
      NSVertex::Vertex<FaceField>, const NSVertex::Vertex<FaceField>*,
      const NSVertex::Vertex<FaceField>&>,
      public MeshXYZConnFAVIterator_base<FaceField>
{
    friend class MeshXYZConnFAVIterator<FaceField>;
	
  public:

    const_MeshXYZConnFAVIterator(FaceField &field_, int cellIndex_,
				 int vertexIndex_)
	: MeshXYZConnFAVIterator_base<FaceField>(field_, cellIndex_,
						 vertexIndex_)
    {
	// empty
    }

    const_MeshXYZConnFAVIterator(const MeshXYZConnFAVIterator<FaceField> &b)
	: MeshXYZConnFAVIterator_base<FaceField>(b)
    {
	// empty
    }

    const_MeshXYZConnFAVIterator &operator++()
    {
	incr(); return *this;
    }
    const_MeshXYZConnFAVIterator operator++(int)
    {
	const_MeshXYZConnFAVIterator tmp = *this; incr(); return tmp;
    }

    typename const_MeshXYZConnFAVIterator::reference operator*()
    { return MeshXYZConnFAVIterator_base<FaceField>::operator*<reference>(); }
    typename const_MeshXYZConnFAVIterator::pointer operator->()
    { return MeshXYZConnFAVIterator_base<FaceField>::operator-><pointer>(); }
};

template<class FaceField>
class NSVertex::Vertex
{
  public:

    typedef typename FaceField::value_type value_type;

  private:
    
    // The iterator for MeshXYZConnFacesAroundVertices
    // uses some private methods of Vertex
    
    friend class MeshXYZConnFAVIterator_base<FaceField>;

    struct FacesAtVertex
    {
	int dat[NUMFACES_PER_VERTX];
	FacesAtVertex(int iv)
	{
	    int vf[NUMVERTICES_PER_CELL][NUMFACES_PER_VERTX] = {
		{0,2,4}, {1,2,4}, {0,3,4}, {1,3,4},
		{0,2,5}, {1,2,5}, {0,3,5}, {1,3,5}
	    };
		
	    Assert(iv >= 0 || iv < NUMVERTICES_PER_CELL);
	    std::copy(vf[iv], vf[iv+1], dat);
	}
    };

    FaceField &field;
    int cellIndex;
    int vertexIndex;

    FacesAtVertex faces;

  private:
    
    class iterator_base;

  public:

    class iterator;
    class const_iterator;

  private:
    
    friend class iterator_base;
    friend class iterator;
    friend class const_iterator;

  private:
    
    class iterator_base
    {
	NSVertex::Vertex<FaceField> &vertex;
	int faceIndex;

      public:

	iterator_base(NSVertex::Vertex<FaceField> &vertex_, int faceIndex_)
	    : vertex(vertex_), faceIndex(faceIndex_)
	{
	    // empty
	};
	
	iterator_base(const iterator_base &fit)
	    : vertex(fit.vertex), faceIndex(fit.faceIndex)
	{
	    // empty
	};

	bool operator==(const iterator_base &rhs) const
	{
	    return &vertex == &rhs.vertex && faceIndex == rhs.faceIndex;
	}
	bool operator!=(const iterator_base &rhs) const
	{
	    return !operator==(rhs);
	}

      protected:

	template<class Reference>
	Reference operator*() const
	{
	    Assert (faceIndex >= 0);
	    Assert (faceIndex < NSVertex::NUMFACES_PER_VERTX);
	    return vertex.field(vertex.cellIndex,
				 vertex.faces.dat[faceIndex]);
	}

	void incr()
	{
	    faceIndex++;
	}

    };

  public:

    class iterator
	: public std::iterator<std::forward_iterator_tag, value_type,
	  typename FaceField::difference_type, value_type*, value_type&>,
	  public iterator_base
    {
      public:

	iterator(NSVertex::Vertex<FaceField> &vertex_, int faceIndex_)
	    : iterator_base(vertex_, faceIndex_)
	{
	    // empty
	}

	iterator(const iterator &it_) : iterator_base(it_) { /* empty */ }
	
	iterator &operator++()
	{
	    incr(); return *this;
	}
	iterator operator++(int)
	{
	    iterator tmp = *this; incr(); return tmp;
	}

	typename iterator::reference operator*() const
	{ return iterator_base::operator*<reference>(); }
    };

    class const_iterator 
	: public std::iterator<std::forward_iterator_tag, value_type,
	  typename FaceField::difference_type, const value_type*,
	  const value_type&>,
	  public iterator_base
    {
      public:
	
	const_iterator(NSVertex::Vertex<FaceField> &vertex_, int faceIndex_)
	    : iterator_base(vertex_, faceIndex_)
	{
	    // empty
	}

	const_iterator(const iterator &it_) : iterator_base(it_) { /* empty */ }
	
	const_iterator &operator++()
	{
	    incr(); return *this;
	}
	const_iterator operator++(int)
	{
	    const_iterator tmp = *this; incr(); return tmp;
	}

	typename const_iterator::reference operator*() const
	{ return iterator_base::operator*<reference>(); }
    };

  public:

    Vertex(FaceField &field_, int cellIndex_, int vertexIndex_)
	: field(field_), cellIndex(cellIndex_), vertexIndex(vertexIndex_),
	 faces(vertexIndex_)
    {
	// empty
    }

    iterator begin()
    {
	return iterator(*this, 0);
    }
    iterator end()
    {
	return iterator(*this, NUMFACES_PER_VERTX);
    }
    
    const_iterator begin() const
    {
	return iterator(const_cast<Vertex&>(*this), 0);
    }
    const_iterator end() const
    {
	return iterator(const_cast<Vertex&>(*this), NUMFACES_PER_VERTX);
    }

  private:

    void assign(int cellIndex_, int vertexIndex_)
    {
	cellIndex = cellIndex_;
	vertexIndex = vertexIndex_;
	faces = FacesAtVertex(vertexIndex);
    }
    
};

} // end namespace rtt_mesh

#endif                          // __mesh_MeshXYZConnFAVIterators_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/MeshXYZConnFAVIterators.hh
//---------------------------------------------------------------------------//
