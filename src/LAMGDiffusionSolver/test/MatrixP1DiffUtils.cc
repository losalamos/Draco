//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/MatrixP1DiffUtils.cc
 * \author Randy M. Roberts
 * \date   Fri Jan 21 10:48:20 2000
 * \brief  Utilities for testing MatrixP1Diff
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MatrixP1DiffUtils.hh"
#include "LAMG/CompressedRowStorage.hh"
#include <algorithm>
#include "ds++/Assert.hh"
#include "c4/global.hh"
#include <iostream>

namespace
{

typedef rtt_LAMG::CompressedRowStorage::IMat1 IMat1;
typedef rtt_LAMG::CompressedRowStorage::DMat1 DMat1;
typedef DMat1::value_type FloatType;
    
void createGlobalCRS(IMat1 &rowPointer_in, IMat1 &colIndex_in, DMat1 &val_in)
{
    // construct the global matrix on process 0

    std::vector<IMat1::value_type> rowPointer(rowPointer_in.begin(),
					      rowPointer_in.end());
    std::vector<IMat1::value_type> colIndex(colIndex_in.begin(),
					    colIndex_in.end());
    std::vector<DMat1::value_type> val(val_in.begin(), val_in.end());
    
    int nentries = colIndex.size();
    Check(val.size() == nentries);
	
    int nrows = rowPointer.size() - 1;

    if (C4::node() == 0)
    {
	// Resize rowPointer removing the last entry.

	rowPointer.resize(nrows);

	for (int source = 1; source < C4::nodes(); source++)
	{
	    int buf[2];
	    C4::Recv(buf, 2, source);
	    int nrows_proc = buf[0];
	    int nentries_proc = buf[1];

	    int *rows_proc = new int[nrows_proc];
	    int *cols_proc = new int[nentries_proc];
	    FloatType *val_proc = new FloatType[nentries_proc];
	    
	    C4::Recv(rows_proc, nrows_proc, source);
	    C4::Recv(cols_proc, nentries_proc, source);
	    C4::Recv(val_proc, nentries_proc, source);

	    // Renumber rows_proc to start at the correct location.

	    for (int i=0; i<nrows_proc; i++)
		rows_proc[i] += nentries;

	    // Resize our arrays.

	    rowPointer.resize(nrows + nrows_proc);
	    colIndex.resize(nentries + nentries_proc);
	    val.resize(nentries + nentries_proc);

	    // Append the proc's data.
	    
	    std::copy(rows_proc, rows_proc + nrows_proc,
		      rowPointer.begin() + nrows);
	    std::copy(cols_proc, cols_proc + nentries_proc,
		      colIndex.begin() + nentries);
	    std::copy(val_proc, val_proc + nentries_proc,
		      val.begin() + nentries);

	    delete [] rows_proc;
	    delete [] cols_proc;
	    delete [] val_proc;

	    nrows += nrows_proc;
	    nentries += nentries_proc;
	}

	// Add the extra rowPointer index to the last entry.

	Ensure(rowPointer.size() == nrows);
	rowPointer.resize(nrows + 1);
	rowPointer[nrows] = nentries;

	Ensure(colIndex.size() == nentries);
	Ensure(val.size() == nentries);

	rowPointer_in.redim(nrows + 1);
	std::copy(rowPointer.begin(), rowPointer.end(), rowPointer_in.begin());

	colIndex_in.redim(nentries);
	std::copy(colIndex.begin(), colIndex.end(), colIndex_in.begin());

	val_in.redim(nentries);
	std::copy(val.begin(), val.end(), val_in.begin());
    }
    else
    {
	int dest = 0;

	int buf[2];

	// Remember the rowPointer vector is nrows + 1 long.
	
	int nrows_proc = rowPointer.size() - 1;
	int nentries_proc = nentries;
	
	buf[0] = nrows_proc;
	buf[1] = nentries_proc;

	C4::Send(buf, 2, dest);

	int *rows_proc = new int[nrows_proc];
	int *cols_proc = new int[nentries_proc];
	FloatType *val_proc = new FloatType[nentries_proc];

	std::copy(rowPointer.begin(), rowPointer.begin() + nrows_proc,
		  rows_proc);
	std::copy(colIndex.begin(), colIndex.end(), cols_proc);
	std::copy(val.begin(), val.end(), val_proc);

	C4::Send(rows_proc, nrows_proc, dest);
	C4::Send(cols_proc, nentries_proc, dest);
	C4::Send(val_proc, nentries_proc, dest);

	delete [] rows_proc;
	delete [] cols_proc;
	delete [] val_proc;
    }
}

void createGlobalRHS(DMat1 &x_in)
{
    // construct the global matrix on process 0
    
    std::vector<DMat1::value_type> x(x_in.begin(), x_in.end());
    
    int nrows = x.size();

    if (C4::node() == 0)
    {
	for (int source = 1; source < C4::nodes(); source++)
	{
	    int buf[1];
	    C4::Recv(buf, 1, source);
	    int nrows_proc = buf[0];

	    FloatType *x_proc = new FloatType[nrows_proc];
	    
	    C4::Recv(x_proc, nrows_proc, source);

	    // Resize our arrays.

	    x.resize(x.size() + nrows_proc);

	    // Append the proc's data.
	    
	    std::copy(x_proc, x_proc + nrows_proc,
		      x.begin() + nrows);

	    delete [] x_proc;

	    nrows += nrows_proc;
	}
	Ensure(x.size() == nrows);
	x_in.redim(nrows);
	std::copy(x.begin(), x.end(), x_in.begin());
    }
    else
    {
	int dest = 0;

	int buf[1];
	int nrows_proc = nrows;
	
	buf[0] = nrows_proc;

	C4::Send(buf, 1, dest);

	FloatType *x_proc = new FloatType[nrows_proc];

	std::copy(x.begin(), x.end(), x_proc);

	C4::Send(x_proc, nrows_proc, dest);

	delete [] x_proc;
    }
}

void doGlobalMultiply(DMat1 &b, const DMat1 &x, const IMat1 &rowPointer,
		      const IMat1 &colIndex, const DMat1 &val)
{
    if (C4::node() == 0)
    {
	int row = 0;
	for (int i = 0; i < colIndex.size(); i++)
	{
	    while (i >= rowPointer[row+1])
		++row;
	    int col = colIndex[i];
	    b[row] += val[i] * x[col];
	    FloatType brow = b[row];
	}
    }
}

void distributeResults(DMat1 &retval, const DMat1 &b)
{
    int nrows = retval.size();
    
    if (C4::node() == 0)
    {
	std::copy(b.begin(), b.begin()+nrows, retval.begin());
	
	for (int source = 1; source < C4::nodes(); source++)
	{
	    int buf[1];
	    C4::Recv(buf, 1, source);
	    int nrows_proc = buf[0];

	    FloatType *b_proc = new FloatType[nrows_proc];
	    
	    // Copy the proc's new data.
	    
	    std::copy(b.begin() + nrows, b.begin() + nrows + nrows_proc,
		      b_proc);

	    C4::Send(b_proc, nrows_proc, source);

	    delete [] b_proc;

	    nrows += nrows_proc;
	}
    }
    else
    {
	int dest = 0;

	int buf[1];
	int nrows_proc = retval.size();
	
	buf[0] = nrows_proc;

	C4::Send(buf, 1, dest);

	FloatType *b_proc = new FloatType[nrows_proc];

	C4::Recv(b_proc, nrows_proc, dest);

	std::copy(b_proc, b_proc + nrows_proc, retval.begin());

	delete [] b_proc;
    }
    
}

} // end unnamed namespace

namespace rtt_LAMGDiffusionSolver_test
{
 
using namespace rtt_LAMGDiffusionSolver;

DMat1 operator*(const MatrixP1Diff &A, const DMat1 &x_in)
{
    // Check for square matrix.

    int nrows_total = A.nrowsTotal();
    Check(A.ncolsTotal() == nrows_total);

    Require(x_in.size() == A.nrows());

    rtt_LAMG::CompressedRowStorage crs = A.crs();

    IMat1 rowPointer = crs.rowPointer();
    IMat1 colIndex = crs.colIndex();
    DMat1 val = crs.val();

    using std::cout;
    using std::endl;

    createGlobalCRS(rowPointer, colIndex, val);

    int nentries_total = crs.colIndex().size();
    C4::gsum(nentries_total);
    
    if (C4::node() == 0)
    {
	Check(rowPointer.size() == nrows_total + 1);
	Check(colIndex.size() == nentries_total);
	Check(val.size() == nentries_total);
    }

    DMat1 x = x_in;

    createGlobalRHS(x);

    if (C4::node() == 0)
    {
	// Check for size of total RHS.
	// This will only be correct on proc 0.
	cout << C4::node() << ": x.size(): " << x.size()
	     << " nrows_total: " << nrows_total << endl;
	Check(x.size() == nrows_total);
    }
    
    DMat1 b(nrows_total, 0.0);

    doGlobalMultiply(b, x, rowPointer, colIndex, val);
    
    int nrows = A.nrows();
    DMat1 retval(nrows);

    distributeResults(retval, b);

    return retval;
}

} // end namespace rtt_LAMGDiffusionSolver_test

//---------------------------------------------------------------------------//
//                              end of MatrixP1DiffUtils.cc
//---------------------------------------------------------------------------//
