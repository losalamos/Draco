//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/PreComputedState.t.hh
 * \author Randy M. Roberts
 * \date   Tue Jan  4 10:13:23 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "PreComputedState.hh"
#include "LAMG/DynamicCompRowStorage.hh"

#include "ds++/Assert.hh"
#include "c4/global.hh"

namespace rtt_LAMGDiffusionSolver
{

template<class MT>
PreComputedState::PreComputedState(const typename MT::FieldConstructor &fCtor,
				   const MT &mesh)
    : nrows_m(mesh.get_ncells())
{
    typename MT::ccif ccif(fCtor);
    diagsToColIndex.resize(ccif.size());

    typename MT::fcdif fcdif(fCtor);
    offDiagsToColIndex.resize(fcdif.size());
    
    // Discern the global offset of row numbers.
    
    calcGlobalRowNumOffset();

    // Determine the Compressed Row Storage Sparsity pattern
    // of the matrix.

    generateCompRowStorage(fCtor, mesh);
}

template<class MT>
void
PreComputedState::generateCompRowStorage(const typename
					 MT::FieldConstructor &fCtor,
					 const MT &mesh)
{
    // Determine the Compressed Row Storage Sparsity pattern
    // of the matrix.
    
    // Create a global numbering of the cells inorder to discern
    // the matrix sparsity pattern.
    // Therefore, obtain a cell-centered field containing the global
    // and local cell numbers.

    typedef typename MT::ccif ccif;
    
    ccif localCellNums(fCtor);

    { // scope
	int cellnum = 0;
	for (ccif::iterator lcit = localCellNums.begin();
	     lcit != localCellNums.end(); ++lcit)
	{
	    *lcit = cellnum; ++cellnum;
	}
    } // scope
    
    ccif globalCellNums(fCtor);

    // Calculate the global cell numbers from the local cell numbers.
    // by adding the global row number offset.

    using std::bind2nd;
    using std::plus;
    
    std::transform(localCellNums.begin(), localCellNums.end(),
		   globalCellNums.begin(), bind2nd(plus<int>(),
						   globalRowNumOffset()));
    Assert(globalCellNums.size() == nrows());

    // We are going to assertain the connectivity (sparse pattern)
    // of the matrix.
    //
    // We know the matrix is represented by a face-centered discontinuous
    // field, so we can...
    //
    // 1) map the global cell numbers on to the faces,
    // 2) do a swap of the faces,
    // 3) and then generate the local part of the sparse matrix pattern.

    // Gather the global and local cell numbers onto the faces.

    typedef typename MT::fcdif fcdif;
    
    fcdif globalCellNumsOnFaces(fCtor);
    MT::gather(globalCellNumsOnFaces, globalCellNums, MT::OpAssign());

    fcdif localCellNumsOnFaces(fCtor);
    MT::gather(localCellNumsOnFaces, localCellNums, MT::OpAssign());

    // Swap the global "unknown" numbers across the faces.
    // This will obtain the column numbers for the sparsity pattern.

    fcdif globalCellNumsAcrossFaces(fCtor);

    // We will use NoElement as a flag to indicate the face
    // on the boundary does not have a matrix element associated
    // with it.
    
    // globalCellNumsAcrossFaces = NoElement;
    MT::swap_faces(globalCellNumsAcrossFaces, globalCellNumsOnFaces);

    typedef typename MT::bsif bsif;
    bsif boundaryFaces(fCtor);
    boundaryFaces = NoElement;
    
    MT::gather(globalCellNumsAcrossFaces, boundaryFaces, MT::OpAssign());

    // The DynamicCompRowStorage will determine the contents of
    // rowPointer and colIndex, and also allow us to create
    // a mapping from the cell-centered diagonal/face-centered
    // off-diagonal matrix representation
    // to the CRS representation.
    
    DynamicCompRowStorage crs;

    addEntries<MT>(crs, localCellNums, globalCellNums,
		   localCellNumsOnFaces, globalCellNumsAcrossFaces);

    setIndexMaps<MT>(crs, localCellNums, globalCellNums,
		     localCellNumsOnFaces, globalCellNumsAcrossFaces);
    
}

    
template<class MT>
void PreComputedState::addEntries(DynamicCompRowStorage &crs,
				  const typename MT::ccif &rowCellNums,
				  const typename MT::ccif &colCellNums,
				  const typename MT::fcdif &rowFaceNums,
				  const typename MT::fcdif &colFaceNums)
{
    // Add diagonal elements to crs.

    crs.add(rowCellNums.begin(), rowCellNums.end(), colCellNums.begin());
    
    // Add off-diagonal elements to crs.

    // The faces on the boundary will contain NoElement
    // as the column number.
    // We do not want to add these faces to the matrix.

    for (typename MT::fcdif::const_iterator rit = rowFaceNums.begin(),
	     cit = colFaceNums.begin(); rit != rowFaceNums.end();
	 ++rit, ++cit)
    {
	if (*cit != NoElement)
	    crs.add(*rit, *cit);
    }

    numEntries_m = crs.numEntries();

    std::cout << C4::node() << ": numEntries_m: " << numEntries_m << std::endl;
}

template<class ITER, class CITER>
void PreComputedState::indices(DynamicCompRowStorage &pat,
			       ITER OutputIter, const ITER OutputEnd,
			       CITER rowIter, CITER colIter) const
{
    // Find the index into columnIndices() for the (local row)/(global column).

    while (OutputIter != OutputEnd)
    {
	if (*colIter == NoElement)
	{
	    *OutputIter = NoElement;
	}
	else
	{
	    *OutputIter = pat.index(*rowIter, *colIter);
	}
	++OutputIter; ++rowIter; ++colIter;
    }
}

template<class MT>
void PreComputedState::setIndexMaps(DynamicCompRowStorage &crs,
				    const typename MT::ccif &rowCellNums,
				    const typename MT::ccif &colCellNums,
				    const typename MT::fcdif &rowFaceNums,
				    const typename MT::fcdif &colFaceNums)
{
    // Set our row pointers to whatever the sparsity pattern says it is.

    { // scope
	
	std::vector<DynamicCompRowStorage::size_type> temp = crs.rowPointers();
    
	std::copy(temp.begin(), temp.end(),
		  std::back_inserter(rowPointer_m));
	std::cout << C4::node()
		  << ": rowPointer_m.size(): " << rowPointer_m.size()
		  << ", nrows() + 1: " << nrows() + 1 << std::endl;
	Assert(rowPointer_m.size() == temp.size());
	Assert(rowPointer_m.size() == nrows() + 1);
    } // end scope
    
    // Set our column indices to whatever the sparsity pattern says it is.
    
    { // scope
	
	std::vector<DynamicCompRowStorage::size_type> temp =
	    crs.columnIndices();
    
	std::copy(temp.begin(), temp.end(),
		  std::back_inserter(colIndex_m));
	Assert(colIndex_m.size() == temp.size());

    } // end scope

    // Fill up the diagsToColIndex with a mapping between
    // cell locations and their associated position within
    // the sparse matrix representation.

    Assert(diagsToColIndex.size() == rowCellNums.size());
    crs.indices(diagsToColIndex.begin(), diagsToColIndex.end(),
		rowCellNums.begin(), colCellNums.begin());

    // Fill up the offDiagsToColIndex with a mapping between
    // face locations and their associated position within
    // the sparse matrix representation.

    Assert(offDiagsToColIndex.size() == rowFaceNums.size());
    indices(crs, offDiagsToColIndex.begin(), offDiagsToColIndex.end(),
	    rowFaceNums.begin(), colFaceNums.begin());
}

} // end namespace rtt_LAMGDiffusionSolver

//---------------------------------------------------------------------------//
//                        end of LAMGDiffusionSolver/PreComputedState.t.hh
//---------------------------------------------------------------------------//
