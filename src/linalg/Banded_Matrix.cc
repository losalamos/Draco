//----------------------------------*-C++-*----------------------------------//
// Banded_Matrix.cc
// Geoffrey M. Furnish
// Thu Mar 12 12:51:51 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "linalg/Banded_Matrix.hh"

#include "c4/SpinLock.hh"

//---------------------------------------------------------------------------//
// Compute and store the info which describes the messaging pattern required
// to perform mathematically correct matrix-vector multiplies.
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::compute_message_routing()
{
// Build up an "occupancy vector" which shows which column positions are
// filled in any row owned by this processor.

    Mat1<int> occupancy( nrt );
    occupancy = 0;

// Now figure out which columns we use.

    for( int n=0; n < N; n++ )
    {
    // Look at diagonal number n.

	for( int c=0; c < nrp; c++ ) {
	    int col = nro + diag_offset[n] + c;
	    if ( 0 <= col && col < nrt )
		occupancy( col ) = 1;
	}
    }

    if (verbose) emit_occupancy_vectors( occupancy );

// Now we need to deduce which processors we will need to send elements from
// our column vector to, as well as the inverse--which processors will be
// sending us data from their column vector.  

// Let's start by getting the occupancy vector for each node communicated to
// each other node.

    Mat2<int> ov( nodes(), nrt );
    {
	Mat1<C4_Req> rreq( nodes() );
	for( int i=0; i < nodes(); i++ ) {
	    if (i == node()) continue;
	    RecvAsync( rreq[i], &ov(i,0), nrt, i );
	}
	for( int i=0; i < nodes(); i++ ) {
	    if (i == node()) continue;
	    Send( &occupancy(0), nrt, i );
	}
    // Now copy ours into place.
	for( int i=0; i < nrt; i++ )
	    ov( node(), i ) = occupancy(i);
    }

// Okay, we have the occupancy vectors for all nodes now, so we should be
// able to start calculating the required messaging patterns.
}

//---------------------------------------------------------------------------//
// For each processor, print out the occupancy vector.  This should look a
// little bit like the fill pattern for the matrix being represented, except
// that diagonals get turned into blocks based on which processors own which
// rows. 
//---------------------------------------------------------------------------//

template<class T, int N>
void Banded_Matrix<T,N>::emit_occupancy_vectors( const Mat1<int>& ov ) const
{
    cout << flush;
    HTSyncSpinLock h;

    for( int i=0; i < nrp; i++ ) {
	cout << node() << "  " << nro+i << ":  ";
	for( int j=0; j < nrt; j++ )
	    cout << ov(j) << "  ";
	cout << endl;
    }
}

//---------------------------------------------------------------------------//
//                              end of Banded_Matrix.cc
//---------------------------------------------------------------------------//
