//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/Cell.hh
 * \author Randy M. Roberts
 * \date   Wed Sep 15 10:12:58 1999
 * \brief  Header file for the Cell class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_Cell_hh__
#define __meshTest_Cell_hh__

#include <set>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

namespace rtt_meshTest
{
 
//===========================================================================//
/*!
 * \class Cell
 *
 * \brief The Cell class is a utility to determine the sizes of
 *        an n-dimensional hexahedron.  Also allows comparison through
 *        operator<().
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Cell 
{

    // NESTED CLASSES AND TYPEDEFS
  public:

    //! Typedef determining the input to the some of the Cell constructors.
    
    typedef std::set<int> VertexList;
    
  private:
    
    // DATA

    int dimension_m;
    std::set<int> vertices_m;
    
  public:

    // CREATORS

    //! Constructor based on VertexList for vertices of a cell.
    //! Throws a std::invalid_argument if not valid.
    
    Cell(int dimension_in, const std::set<int> &vertices_in)
	: dimension_m(dimension_in), vertices_m(vertices_in)
    {
	checkValid();
    }

    //! Constructor based on vector of ints for vertices of a cell.
    //! Throws a std::invalid_argument if not valid.
    
    Cell(int dimension_in, const std::vector<int> &vertices_in)
	: dimension_m(dimension_in), vertices_m(vertices_in.begin(),
						vertices_in.end())
    {
	checkValid();
    }

    //! Constructor based on iterators over vertices of a cell.
    //! Throws a std::invalid_argument if not valid.
    
    template<class InputIterator>
    Cell(int dimension_in, InputIterator first, InputIterator last)
	: dimension_m(dimension_in), vertices_m(first, last)
    {
	checkValid();
    }

    //! Destructor
    
    ~Cell() { /* empty */ }

    // MANIPULATORS
    
    // ACCESSORS

    //! Check to make sure that the vertices given to the cell constructor
    //! are valid.  Throws a std::invalid_argument if not valid.
    
    void checkValid() const
    {
	if (vertices().size() != numVertices())
	{
	    std::ostringstream ost;
	    ost << "Cell constructor: size of vertex set"
		<< "(" << vertices_m.size() << ")"
		<< " != "
		<< "correct number of vertices"
		<< "(" << numVertices() << ")"
		<< " for " << dimension() << " dimensions";
	    throw std::invalid_argument(ost.str());
	}
    }

    //! Returns the number of vertices per cell.
    
    int numVertices() const
    {
	// 2 ^ Dim
	return 1 << dimension();
    }

    //! Returns the number of faces per cell.
    
    int numFaces() const
    {
	// 2 * Dim
	return 2 * dimension();
    }

    //! Returns the number of edges per cell.
    
    int numEdges() const
    {
	// Dim * 2 ^ (Dim -1)
	return dimension() * (1 << (dimension()-1));
    }
    
    //! Returns the number of faces surrounding a vertex of the cell.
    
    int numFacesArroundVertex() const
    {
	return dimension();
    }

    //! Returns the dimension of the cell.

    int dimension() const { return dimension_m; }

    //! Returns the list of vertices in this cell.
    
    const VertexList &vertices() const { return vertices_m; }

    //! Determines an ordering between cells.
    //! Throws std::length_error if the cells are of different dimensions.
    
    bool operator<(const Cell &c) const
    {
	if (dimension() != c.dimension())
	    throw std::length_error(
		"Comaring two cells of different dimensions.");

	// Lexicographic comparison of sets of ints.
	return vertices() < c.vertices();
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_Cell_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/Cell.hh
//---------------------------------------------------------------------------//
