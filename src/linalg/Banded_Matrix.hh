//----------------------------------*-C++-*----------------------------------//
// Banded_Matrix.hh
// Geoffrey M. Furnish
// Thu Mar 12 12:51:51 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __linalg_Banded_Matrix_hh__
#define __linalg_Banded_Matrix_hh__

#include "ds++/Mat.hh"
using namespace dsxx;

//===========================================================================//
// class Banded_Matrix - 

// 
//===========================================================================//

template<class T, int N>
class Banded_Matrix
{
    int nro;			// row offset.
    int nrp, nrt;		// # rows this processor, total.

    Mat2<T> data;
    int diag_offset[N];

    bool verbose;

  public:
    template<class Decomposition_DB>
    Banded_Matrix( const Decomposition_DB& ddb, int *doff )
	: nro( ddb.row_offset() ),
	  nrp( ddb.nrows_this_processor() ),
	  nrt( ddb.nrows_total() ),

	  data( nrp, N ),

	  verbose( ddb.verbose() )
    {
    // Store away the diagonal offsets.
	for( int i=0; i < N; i++ )
	    diag_offset[i] = doff[i];

	compute_message_routing();
    }

    void compute_message_routing();

    template<class ColumnVector>
    void multiply( const ColumnVector& x, ColumnVector& b );

// Routines for data access. ...

    T& operator()( int row, int n ) { return data(row,n); }
    const T& operator()( int row, int n ) const { return data(row,n); }

    void emit_occupancy_vectors( const Mat1<int>& ov ) const;
};

#endif                          // __linalg_Banded_Matrix_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/Banded_Matrix.hh
//---------------------------------------------------------------------------//
