//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/DynamicCompRowStorage.hh
 * \author Randy M. Roberts
 * \date   Wed Jan  5 15:41:19 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMG_DynamicCompRowStorage_hh__
#define __LAMG_DynamicCompRowStorage_hh__

#include <vector>
#include <map>
#include <limits>
#include <algobase>

namespace rtt_LAMG
{
 
//===========================================================================//
/*!
 * \class DynamicCompRowStorage
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class DynamicCompRowStorage 
{

    // NESTED CLASSES AND TYPEDEFS

  public:
    
    typedef long int size_type;

  private:
    
    typedef std::map<size_type, size_type> Columns;
    typedef std::map<size_type, Columns> Rows;
    
    // DATA


    // A data structure consisting of the actual sparsity pattern.
    // Mutable so we can call renumber() from const methods.
    
    mutable Rows rows;

    // A flag that detemines if the data structure needs renumbering.
    
    mutable bool dirty;

    size_type  numEntries_m;
    
    // The extrema of the rows and columns stored in pattern.
    
    size_type  rowMax;
    size_type  colMax;
    size_type  rowMin;
    size_type  colMin;

  public:

    // CREATORS

    DynamicCompRowStorage()
	: numEntries_m(0), dirty(false),
	  rowMax(std::numeric_limits<size_type>::min()),
	  colMax(std::numeric_limits<size_type>::min()),
	  rowMin(std::numeric_limits<size_type>::max()),
	  colMin(std::numeric_limits<size_type>::max())
    {
	// empty
    }

    // MANIPULATORS
    
    void add(size_type row, size_type col);

    template<class ITER>
    inline void add(ITER rowIter, const ITER rowEnd, ITER colIter);
    
    // ACCESSORS

    size_type numRows() const { return rowMax + 1; }
    size_type numCols() const { return colMax + 1; }
    size_type numEntries() const { return numEntries_m; }
    
    // Returns the vector of locations into columnIndices() vector
    // that refer to the given row.
    
    std::vector<size_type> rowPointers() const;

    // Returns a vector, of size numeEntries, of column numbers.
    
    std::vector<size_type> columnIndices() const;

    // Find the index into columnIndices() for the (local row)/(global column).

    size_type index(size_type row, size_type col) const;

    // Find the index into columnIndices() for the (local row)/(global column).

    template<class ITER, class CITER>
    inline void indices(ITER OutputIter, const ITER OutputEnd,
			CITER rowIter, CITER colIter) const;

  private:
    
    // IMPLEMENTATION

    void setExtrema(size_type row, size_type col)
    {
	rowMax = std::max(row, rowMax);
	colMax = std::max(col, colMax);
	rowMin = std::min(row, rowMin);
	colMin = std::min(col, colMin);
    }

    void renumber() const;
    
};

template<class ITER>
void DynamicCompRowStorage::add(ITER rowIter, const ITER rowEnd, ITER colIter)
{
    // Add the (local row)/(global column) pair to the
    // sparsity pattern generator.
	
    while (rowIter != rowEnd)
    {
	add(*rowIter, *colIter);
	++rowIter; ++colIter;
    }
}

template<class ITER, class CITER>
void DynamicCompRowStorage::indices(ITER OutputIter, const ITER OutputEnd,
				    CITER rowIter, CITER colIter) const
{
    // Find the index into columnIndices() for the (local row)/(global column).

    while (OutputIter != OutputEnd)
    {
	*OutputIter = index(*rowIter, *colIter);
	++OutputIter; ++rowIter; ++colIter;
    }
}

} // end namespace rtt_LAMG

#endif // __LAMG_DynamicCompRowStorage_hh__

//---------------------------------------------------------------------------//
// end of LAMG/DynamicCompRowStorage.hh
//---------------------------------------------------------------------------//
