//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/PreComputedState.hh
 * \author Randy M. Roberts
 * \date   Tue Jan  4 10:13:23 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_PreComputedState_hh__
#define __LAMGDiffusionSolver_PreComputedState_hh__

#include "ds++/Assert.hh"

#include <vector>
#include <limits>

// Forward Reference

namespace rtt_LAMG
{
class DynamicCompRowStorage;
} // end namespace rtt_LAMG

namespace rtt_LAMGDiffusionSolver
{

// Forward Reference

class PCS_iterator;
 
//===========================================================================//
/*!
 * \class PreComputedState
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class PreComputedState 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef rtt_LAMG::DynamicCompRowStorage DynamicCompRowStorage;
    
    enum {NoElement=-999};

  public:

    class matrix_iterator;
    typedef matrix_iterator const_matrix_iterator;

  private:

    // DATA

    // number of rows on this processor
    
    int nrows_m;
    
    int globalRowNumOffset_m;

    int numEntries_m;
    
    // cellsToColIndex allows easy mapping of the diagonal part
    // of the matrix from cell-centerd to CRS (compressed row storage)
    
    std::vector<int>  diagsToColIndex;
    
    // facessToColIndex allows easy mapping of the off-diagonal part
    // of the matrix from face-centerd to CRS (compressed row storage)
    
    std::vector<int> offDiagsToColIndex;
    
    std::vector<int> rowPointer_m;
    std::vector<int> colIndex_m;
    
  public:

    // CREATORS

    template<class MT>
    PreComputedState(const typename MT::FieldConstructor &fCtor_in,
		     const MT &mesh);

    // MANIPULATORS
    
    matrix_iterator matrix_begin();
    matrix_iterator matrix_end();


    // ACCESSORS

    int nrows() const { return nrows_m; }
    int numEntries() const { return numEntries_m; }
    int globalRowNumOffset() const { return globalRowNumOffset_m; }
    
    const std::vector<int> &rowPointer() const { return rowPointer_m; }
    const std::vector<int> &colIndex() const { return colIndex_m; }

    template<class ITER1, class ITER2>
    inline std::vector<double> val(ITER1 diagonalIter,
				   ITER1 diagonalEnd,
				   ITER2 offDiagIter,
				   ITER2 offDiagEnd) const;

    const_matrix_iterator matrix_begin() const;
    const_matrix_iterator matrix_end() const;

  private:
    
    // IMPLEMENTATION

    template<class MT>
    void generateCompRowStorage(const typename MT::FieldConstructor &fCtor_in,
				const MT &mesh);


    template<class ITER, class CITER>
    void indices(DynamicCompRowStorage &pat,
		 ITER OutputIter, const ITER OutputEnd,
		 CITER rowIter, CITER colIter) const;

    void calcGlobalRowNumOffset();

    template<class MT>
    void addEntries(DynamicCompRowStorage &crs,
		    const typename MT::ccif &rowCellNums,
		    const typename MT::ccif &colCellNums,
		    const typename MT::fcdif &rowFaceNums,
		    const typename MT::fcdif &colFaceNums);

    template<class MT>
    void setIndexMaps(DynamicCompRowStorage &crs,
		      const typename MT::ccif &rowCellNums,
		      const typename MT::ccif &colCellNums,
		      const typename MT::fcdif &rowFaceNums,
		      const typename MT::fcdif &colFaceNums);

  public:

    // Nested class definition.
    
    class matrix_iterator 
    {

	// NESTED CLASSES AND TYPEDEFS

	const PreComputedState &pcs;
	int entry;
	int row_m;
	int col_m;

	friend class PreComputedState;
	
	// DATA
    
	// CREATORS

      private:

	// All creator are private.
	
	matrix_iterator(const PreComputedState &pcs_in, int entry_in)
	    : pcs(pcs_in), entry(entry_in)
	{
	    Require(entry_in >= 0 && entry_in <= pcs.numEntries());
	    
	    row_m = 0;
	    if (entry_in < pcs.numEntries())
		calcRowColumn();

	    Ensure(entry >= 0 && entry <= pcs.numEntries());
	}

      public:
    
	// MANIPULATORS

	matrix_iterator &operator++()
	{
	    Check(entry >= 0 && entry <= pcs.numEntries());
	    
	    if (entry < pcs.numEntries())
	    {
		++entry;
		calcRowColumn();
	    }
	    else
	    {
		entry = pcs.numEntries();
	    }
	    return *this;
	}
    
	matrix_iterator operator++(int)
	{
	    Check(entry >= 0 && entry <= pcs.numEntries());

	    matrix_iterator tmp = *this;
	    ++*this;
	    return tmp;
	}
    
	// ACCESSORS

	bool operator==(const matrix_iterator &rhs) const
	{
	    Check(entry >= 0 && entry <= pcs.numEntries());
	    Check(rhs.entry >= 0 && rhs.entry <= rhs.pcs.numEntries());
	    
	    return &pcs == &rhs.pcs && entry == rhs.entry;
	}
    
	bool operator!=(const matrix_iterator &rhs) const
	{
	    return !(*this == rhs);
	}
    
	int row() const { return row_m; }
	int col() const { return col_m; }

      private:
    
	// IMPLEMENTATION

	void calcRowColumn()
	{
	    while (entry >= pcs.rowPointer()[row_m+1])
		row_m++;
	    col_m = pcs.colIndex()[entry];
	}
    };
};

template<class ITER1, class ITER2>
std::vector<double> PreComputedState::val(ITER1 diagonalIter,
					  ITER1 diagonalEnd,
					  ITER2 offDiagIter,
					  ITER2 offDiagEnd) const
{
    // Retrieve the matrix elements stored in the compressed row storage
    // format.

    std::vector<double> retval(numEntries(),
			       std::numeric_limits<double>::min());

    Assert(std::distance(diagonalIter, diagonalEnd) == diagsToColIndex.size());

    { // scoping

	std::vector<int>::const_iterator dciit = diagsToColIndex.begin();

	while (diagonalIter != diagonalEnd)
	{
	    Assert(*dciit >= 0);
	    Assert(*dciit < numEntries());
	    retval[*dciit] = *diagonalIter; ++dciit; ++diagonalIter;
	}

    } // end of scoping

    Assert(std::distance(offDiagIter, offDiagEnd) == offDiagsToColIndex.size());

    { // scoping
	
	std::vector<int>::const_iterator odciit = offDiagsToColIndex.begin();

	while (offDiagIter != offDiagEnd)
	{
	    // NoElement is a special symbol for no matrix element.
	    // This will happen at the faces along the boundary.
	    
	    if (*odciit != NoElement)
	    {
		Assert(*odciit >= 0);
		Assert(*odciit < numEntries());
		retval[*odciit] = *offDiagIter;
	    }
	    ++odciit; ++offDiagIter;
	}

    } // end of scoping

    return retval;
}

} // end namespace rtt_LAMGDiffusionSolver

#endif                          // __LAMGDiffusionSolver_PreComputedState_hh__

//---------------------------------------------------------------------------//
//                              end of LAMGDiffusionSolver/PreComputedState.hh
//---------------------------------------------------------------------------//
