//----------------------------------*-C++-*----------------------------------//
// T3_MatVec.cc
// Dave Nystrom
// Thu Oct  2 14:21:13 MDT 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "t3/T3_MatVec.hh"

#include "t3/Test_Prob.hh"

#include "c4/global.hh"
#include "c4/C4_Req.hh"
#include "c4/SpinLock.hh"
using namespace C4;

#include <stdio.h>

//---------------------------------------------------------------------------//
// Destructor.
//---------------------------------------------------------------------------//

template<class T>
T3_MatVec<T>::~T3_MatVec() {}

//---------------------------------------------------------------------------//
// Constructor.
//---------------------------------------------------------------------------//

template<class T>
T3_MatVec<T>::T3_MatVec( Test_Prob *p )
    : PCG_MatVec<T>(), prob(p),
      ncps(nodes), goffs(nodes),
      its(0)
{
    ncps = 0;
    goffs = 0;
    ncps[node] = prob->ncp;
    goffs[node] = prob->goff;

//     char buf[80];

    for( int n=0; n < nodes; n++ ) {
	C4::gsum( ncps[n] );
	C4::gsum( goffs[n] );
    }

//     {
// 	HTSyncSpinLock h;

// 	for( int n=0; n < nodes; n++ ) {
// 	    sprintf( buf, "node %d, ncps[%d]=%d goffs[%d]=%d",
// 		     node, n, ncps[n], n, goffs[n] );
// 	    cout << buf << endl;
// 	}
//     }
}

//---------------------------------------------------------------------------//
// Evaluate matrix-vector product.
//---------------------------------------------------------------------------//

template<class T>
void T3_MatVec<T>::MatVec( Mat1<T>& b, Mat1<T>&x )
{
    its++;
    char buf[80];

    int ncp = prob->ncp;
    int nct = prob->nct;
    Mat2<double>& A = prob->A;

    {
	Mat1<double> xx( nct );
	Mat1<C4_Req> /*sreq(nodes-1),*/ rreq(nodes-1);

	Assert( x.nx() == ncp );

    // Post the async recvs.

	for( int m=0, n=0; n < nodes; n++ ) {
	    if (n == node) continue;
	    RecvAsync( rreq[m++], &xx( goffs[n] ), ncps[n]*sizeof(double),
		       n, C4_double_ptr_Tag );
	}

    // Post the async sends.

	for( int m=0, n=0; n < nodes; n++ ) {
	    if (n == node) continue;
	//	    SendAsync( sreq[m++], &x( 0 ), ncp*sizeof(double), n, 1
	//);
	    Send( &x(0), ncp, n );
	}

    // Copy our own data into place.

	for( int i=0; i < ncp; i++ )
	    xx[ goffs[node] + i ] = x[i];

    // Wait on the recv's.

	for( int m=0, n=0; n < nodes; n++ ) {
	    if (n == node) continue;
	    rreq[m++].wait();
	}

// 	for( int i=0; i < nct; i++ ) {
// 	    sprintf( buf, "node %d xx[%d]=%lf", node, i, xx[i] );
// 	    cout << buf << endl;
// 	}

    // Perform multiplication.

	for( int i=0; i < ncp; i++ ) {
	    b[i] = 0.;
	    for( int j=0; j < nct; j++ )
		b[i] += A(i,j) * xx(j);
	}

    }
}

//---------------------------------------------------------------------------//
//                              end of T3_MatVec.cc
//---------------------------------------------------------------------------//
