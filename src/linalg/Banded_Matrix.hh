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

#include "c4/C4_Req.hh"

//===========================================================================//
// class Banded_Matrix - 

// 
//===========================================================================//

template<class T, int N>
class Banded_Matrix
{
    int nro;			// row offset.
    int nrp, nrt;		// # rows this processor, total.

    dsxx::Mat2<T> data;		// The actual non-zero elements in the matrix.
    int diag_offset[N];		// Band offset from diagonal.

// For a given row, indicate the range of active bands

    dsxx::Mat1<int> start_band, nbands;

// Receive buffer for the data from the column vector which will be sent to
// us during the matrix vector multiply routine.

    dsxx::Mat1<T> cvdata;

// Mapping from band, row indices to the index in cvdata which holds the
// corresponding element from the vector.

    dsxx::Mat2<int> cvindex;

    bool verbose;		// Be chatty?

    int receivers;
    dsxx::Mat1<int> receiver_node, receiver_nels;
    dsxx::Mat2<int> receiver_ndxs;
    dsxx::Mat2<T>   receiver_data;
    dsxx::Mat1<C4::C4_Req> sreq;

    int self_receiver_id, self_sender_id;

    int senders;
    dsxx::Mat1<int> sender_node, sender_nels, sender_head;
    dsxx::Mat1<C4::C4_Req> rreq;

    struct rcvr_tag {
	int node;
	int nels;
	dsxx::Mat1<int> indexes;

	rcvr_tag( int node_, int nels_, int *indexes_ )
	    : node(node_), nels(nels_), indexes(nels)
	{
	    std::copy( indexes_, indexes_ + nels, indexes.begin() );
	}
    };

  public:
    template<class Decomposition_DB>
    Banded_Matrix( const Decomposition_DB& ddb, const int *doff )
	: nro( ddb.row_offset() ),
	  nrp( ddb.nrows_this_processor() ),
	  nrt( ddb.nrows_total() ),

	  data( N, nrp ),

	  start_band(nrp), nbands(nrp),

	  cvindex( N, nrp ),

	  verbose( ddb.verbose() ),

	  self_receiver_id(-1), self_sender_id(-1)
    {
    // Store away the diagonal offsets.
	for( int i=0; i < N; i++ )
	    diag_offset[i] = doff[i];

	compute_message_routing();
    }

    void compute_message_routing();

    Banded_Matrix& operator=( T x ) { data = x; return *this; }

    void initiate_sends( const T *pd );
    void complete_sends();

    void initiate_receives();
    void complete_receives();

    template<class ColumnVector>
    void multiply( const ColumnVector& x, ColumnVector& b );

// Routines for data access. ...

    T& operator()( int row, int n ) { return data(n,row); }
    const T& operator()( int row, int n ) const { return data(n,row); }

    void emit_occupancy_vectors( const dsxx::Mat1<int>& ov ) const;
    void emit_receiver_structs() const;
    void emit_sender_structs() const;
};

#endif                          // __linalg_Banded_Matrix_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/Banded_Matrix.hh
//---------------------------------------------------------------------------//
