//----------------------------------*-C++-*----------------------------------//
// dump.cc
// Maurice LeBrun
// Thu Feb  9 15:53:47 1995
//---------------------------------------------------------------------------//
// @> Functions for dumping various things in various ways.
//---------------------------------------------------------------------------//

#include "dump.hh"
#include "c4/SpinLock.hh"
using namespace C4;

template<class T>
void print_Mat( const Mat1<T>& a, String name )
{
    int node = C4_node();

    NODE0 cout << "Printing Mat1<T> " << name << " :\n";
    {
	HTSyncSpinLock h;

	cout << "node : " << node << endl;
	int s = a.nx();
	cout << " size=" << s << endl;
	for( int i=0; i < s; i++ )
	    cout << a[i] << ' ';
	cout << endl;
    }
}

template<class T>
void print_Mat( const Mat2<T>& m, String name )
{
    int node = C4::node();
//     int nfu, ncv;
//     m.elements( nfu, ncv );

    int nx = m.nx();
    int ny = m.ny();

    NODE0 cout << "Printing Mat2<T> " << name << " :\n";

    for( int r=0; r < nx; r++ ) {

    //	NODE0 cout << "j=" << j << endl;

	{
	    HTSyncSpinLock h;

	    for( int c=0; c < ny; c++ ) {
		cout.width(6);
		cout << m(r,c) << ' ';
	    }
	}

	NODE0 cout << endl << flush;
    }
}

template<class T>
void print2_Mat( const Mat2<T>& m, String name )
{
// Because of the sucky buffering properties on the Origin 2000, we're going
// to have to got to special effort here.

    cout << flush;
    C4::gsync();

// Even that won't be guaranteed to flush the output from all the procs
// before we start printing below.  Sheesh.

    int node = C4::node();
    int nodes = C4::nodes();

    if (node == 0) {

	cout << "Printing a Mat2(" << m.nx()*nodes << ","
	     << m.ny() << "), " << name << endl;

    // Print my stuff.
	int row=0;
	for( int r=0; r < m.nx(); r++, row++ )
	{
	    cout.width(4);
	    cout << row << ": ";
	    for( int c=0; c < m.ny(); c++ ) {
		cout.width(6);
		cout << m(r,c) << ' ';
	    }
	    cout << endl;
	}

    // Now cycle through the other processors, fetching their stuff and then
    // printing it.

	Mat2<T> mm( m.nx(), m.ny() );

	for( int n=1; n < nodes; n++ )
	{
	    Send( 1, n );
	    Recv( &mm(0,0), m.nx()*m.ny(), n );

	    for( int r=0; r < m.nx(); r++, row++ )
	    {
		cout.width(4);
		cout << row << ": ";
		for( int c=0; c < m.ny(); c++ ) {
		    cout.width(6);
		    cout << mm(r,c) << ' ';
		}
		cout << endl;
	    }
	}
    }
    else {
    // Wait on node 0 to tell me he's ready.
	int dummy;
	Recv( dummy, 0 );

    // Send him my data.
	Send( &m(0,0), m.nx()*m.ny(), 0 );
    }

    C4::gsync();
}

//---------------------------------------------------------------------------//
// This function dumps a distributed 3-d matrix out to an amorphous data file.
// Floowing the example of the following function, this one gathers data to
// one node for bulk output.  Still loop over the third index though, so that
// the memory requirements on the master node don't cause it to fail.
//---------------------------------------------------------------------------//

template<class T>
void ADdump_Mat( ADFile *f, Mat3<T>& m, String id )
{
    int ncu, nfu, ncv, ncw;
    m.elements(nfu,ncv,ncw);
    ADKey adk;
    int node = C4::node(), nodes = C4::nodes(), group = C4::group();
    ncu = nfu * nodes;

    int i, j;
    Mat2<T> mm( ncu, ncv );

    if (!node) {
	sprintf( adk.s, "Mat3f(%d,%d,%d) %s",
		 ncu, ncv, ncw, &id[0] );
	f->Start_new_entry(adk);
    }

    for( int k=0; k < ncw; k++ ) {

    // Now gather in this plane.

	if (!node) {
	    for( j=0; j < ncv; j++ )
		for( i=0; i < nfu; i++ )
		    mm(i,j) = m(i,j,k);

	    if (nodes > 1) {

		Mat2<T> buf( nfu, ncv );

	    // Now get the data from the other nodes and stuff it into the
	    // staging area also.

		for( int n=1; n < nodes; n++ ) {

		    throw "Bogus";
// 		    C4_Recv( &buf(0,0), sizeof(T)*nfu*ncv,
// 			     n, 200, group );

		    for( j=0; j < ncv; j++ )
			for( i=0; i < nfu; i++ )
			    mm(nfu*n+i,j) = buf(i,j);
		}
	    }

	// Now write it out and finish up.

	    f->write( (void *) &mm(0,0), sizeof(T)*ncu*ncv );
	}
	else
	    throw "Bogus";
// 	    C4_Send( &m(0,0,k), sizeof(T) * nfu * ncv,
// 		     0, 200, group );

	gsync();
    }

    if (!node)
	f->End_of_entry();
}

//---------------------------------------------------------------------------//
// Dump a distributed 2-d matrix out to the amorphous data file.  This one
// operates by collecting the data to one node (0) and writing it out in one
// shot, which should be the fastest way to do it.
//---------------------------------------------------------------------------//

template<class T>
void ADdump_Mat( ADFile *f, Mat2<T>& m, String id )
{
    int ncu, nfu, ncv;
    m.elements( nfu, ncv );

    int i, j;
    int node = C4::node(), nodes = C4::nodes(), group = C4::group();

    ncu = nfu * nodes;

    Mat2<T> mm(ncu,ncv);

    if (!node) {
	ADKey adk;
	sprintf( adk.s, "Mat2f(%d,%d) %s",
		 ncu, ncv, &id[0] );
	f->Start_new_entry( adk );
    }

    if (!node) {

    // Copy our data to the staging area.

	for( j=0; j < ncv; j++ )
	    for( i=0; i < nfu; i++ )
		mm(i,j) = m(i,j);

	if (nodes > 1) {

	    Mat2<T> buf( nfu, ncv );

	// Now get the data from the other nodes and stuff it into the staging
	// area also.

	    for( int n=1; n < nodes; n++ ) {

		throw "Bogus";
// 		C4_Recv( &buf(0,0), sizeof(T) * nfu * ncv,
// 			 n, 250, group );

		for( j=0; j < ncv; j++ )
		    for( i=0; i < nfu; i++ )
			mm(nfu*n+i,j) = buf(i,j);
	    }
	}

    // Now write it out and finish up.

	f->write( (void *) &mm(0,0), sizeof(T)*ncu*ncv );
	f->End_of_entry();
    }
    else
	throw "Bogus";
// 	C4_Send( &m(0,0), sizeof(T) * nfu * ncv, 0, 250, group );

    gsync();
}

//---------------------------------------------------------------------------//
//                              end of dump.cc
//---------------------------------------------------------------------------//
