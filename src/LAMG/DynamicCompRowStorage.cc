//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/DynamicCompRowStorage.cc
 * \author Randy M. Roberts
 * \date   Wed Jan  5 15:41:19 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DynamicCompRowStorage.hh"

#include <iostream>

#if 1

#include "ds++/Assert.hh"
#include "c4/global.hh"

#else

#include <sstream>
#include <stdexcept>

#define Assert(arg) assert(arg, __LINE__)

namespace {

void assert(bool assertion, int line)
{
    std::ostringstream strline;
    strline << line;
    if (!assertion)
	throw std::runtime_error("caught exception on line " + strline.str());
}

}

#endif

namespace rtt_LAMG
{

void DynamicCompRowStorage::add(size_type row, size_type col)
{
    Assert(row >= 0);
    Assert(col >= 0);
    
    using std::pair;
    
    setExtrema(row, col);

    // Try to insert.
    // If it already exists, that's OK.  It will act just like
    // a find() instead.
    
    pair<Rows::iterator, bool> rowSuccess =
	rows.insert(Rows::value_type(row, Columns()));
    Columns &columnsForRow = rowSuccess.first->second;

    pair<Columns::iterator, bool> colSuccess
	= columnsForRow.insert(Columns::value_type(col,0));

    // See if we actually inserted anything new.
    
    if (colSuccess.second == true)
    {
	++numEntries_m;

	// We need to renumber the data structure.

	dirty = true;
    }
}

DynamicCompRowStorage::size_type
DynamicCompRowStorage::index(size_type row, size_type col) const
{
    // Find the index into columnIndices() for the (local row)/(global column).

    Assert(row >= rowMin);
    Assert(row <= rowMax);
    Assert(col >= colMin);
    Assert(col <= colMax);

    // Renumber if necessary

    renumber();
    
    Rows::const_iterator rit = rows.find(row);
    Assert(rit != rows.end());

    // The map of columns to indices into columnIndex
    // is the mapped_type (i.e. second of the pair)
    
    const Columns &columnsForRow = rit->second;

    Columns::const_iterator cit = columnsForRow.find(col);

    Assert(cit != columnsForRow.end());
    
    // The index into columnIndex
    // is the mapped_type (i.e. second of the pair)

    int colInd = cit->second;

    Assert(colInd >= 0);
    Assert(colInd < numEntries_m);
	
    return colInd;
}

std::vector<DynamicCompRowStorage::size_type>
DynamicCompRowStorage::rowPointers() const
{
    // Returns the vector of locations into columnIndices() vector
    // that refer to the given row.
    
    // The container must be one larger than the number of rows;
    // therefore, the number of rows, max row number + 1, +1,
    // or max row number + 2.  (Remember we start at row zero.)
    
    std::vector<size_type> retval(rowMax+2);

    int lastRow = -1;
    int lastRowPtr = 0;
    for (Rows::const_iterator rit = rows.begin();
	 rit != rows.end(); ++rit)
    {
	// The row is the key, i.e. first of the pair.
	
	int row = rit->first;

	Assert(row > lastRow);
	Assert(row <= rowMax);
	for (int ir = lastRow+1; ir <= row; ++ir)
	    retval[ir] = lastRowPtr;
	
	// The map of columns to indices into columnIndex
	// is the mapped_type (i.e. second of the pair)
	
	const Columns &columnsForRow = rit->second;
	lastRowPtr += columnsForRow.size();
	lastRow = row;
    }

    Assert(lastRowPtr == numEntries_m);
    retval[rowMax+1] = numEntries_m;

    return retval;
}

void DynamicCompRowStorage::renumber() const
{
    if (!dirty)
	return;

    int index = 0;
    for (Rows::iterator rit = rows.begin(); rit != rows.end(); ++rit)
    {
	// The map of columns to indices into columnIndex
	// is the mapped_type (i.e. second of the pair)
    
	Columns &columnsForRow = rit->second;
	for (Columns::iterator cit = columnsForRow.begin();
	     cit != columnsForRow.end(); ++cit)
	{
	    // The index into columnIndex
	    // is the mapped_type (i.e. second of the pair)

	    cit->second = index++;
	}
    }
    Assert(index == numEntries_m);

    dirty = false;
}

std::vector<DynamicCompRowStorage::size_type>
DynamicCompRowStorage::columnIndices() const
{
    // Returns a vector, of size numeEntries, of column numbers.
    
    std::vector<size_type> retval(numEntries_m);

    int index = 0;
    for (Rows::const_iterator rit = rows.begin(); rit != rows.end(); ++rit)
    {
	// The map of columns to indices into columnIndex
	// is the mapped_type (i.e. second of the pair)
    
	const Columns &columnsForRow = rit->second;
	for (Columns::const_iterator cit = columnsForRow.begin();
	     cit != columnsForRow.end(); ++cit)
	{
	    // The column is the key, i.e. first of the pair.
	    retval[index++] = cit->first;
	}
    }
    return retval;
}

} // end namespace rtt_LAMG

//---------------------------------------------------------------------------//
//                              end of DynamicCompRowStorage.cc
//---------------------------------------------------------------------------//
