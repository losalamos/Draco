//----------------------------------*-C++-*----------------------------------//
// global_nx.cc
// Maurice LeBrun
// Wed Feb  1 16:43:11 1995
//---------------------------------------------------------------------------//
// @> Global C4 functions for the paragon
//---------------------------------------------------------------------------//

#include "DynArray.hh"
#include "Assert.hh"

//---------------------------------------------------------------------------//
// Miscellaneous

void C4_Init( int& argc, char **& argv ) {}

void C4_Finalize()
{
// Get all those annoying termination messages to collect at the bottom.

    gsync();
}

int C4_node()
{
    static int node = -1;
    if (node == -1)
	node = (int) mynode();
    return node;
}

int C4_nodes()
{
    return (int) numnodes();
}

int C4_group()
{
    return (int) myptype();
}

void C4_gsync()
{
    gsync();
}

//---------------------------------------------------------------------------//
// NX send/receive calls
//
// Synchronous:
//    void csend( long type, char *buf, long count, long node, long ptype );
//    void crecv( long typesel, char *buf, long count );
//
// Asynchronous:
//    long isend( long type, char *buf, long count, long node, long ptype );
//    long irecv( long typesel, char *buf, long count );
//	(both return message id)
//    msgdone(mid);
//    msgwait(mid);
//    msgignore(mid);
//
// Extended receive:
//    void crecvx( long typesel, char *buf, long count, long nodesel, long
//		   ptypesel, long info[]);
//    long irecvx( long typesel, char *buf, long count, long nodesel, long
//		   ptypesel, long info[]);
//	(returns message id)

//---------------------------------------------------------------------------//
// Perform a normal (blocking) send.
//---------------------------------------------------------------------------//

int C4_Send( void *buf, int size, int dest, int tag, int group )
{
    csend( tag, (char *) buf, size, dest, group );
    return C4_SUCCESS;
}

//---------------------------------------------------------------------------//
// Perform a normal (blocking) receive.
//---------------------------------------------------------------------------//

int C4_Recv( void *buf, int size, int source, int tag, int group )
{
    if (source == C4_Any_Source)
	crecv( tag, (char *) buf, size );
    else
	crecvx( tag, (char *) buf, size, source, group, msginfo );

    return infocount();
}

//---------------------------------------------------------------------------//
// Perform a non blocking send.
//---------------------------------------------------------------------------//

C4_Req C4_SendAsync( void *buf, int size, int dest, int tag, int group )
{
    C4_Req r;
    r.mid = isend( tag, (char *) buf, size, dest, group );
    r.set();
    return r;
}

//---------------------------------------------------------------------------//
// Perform a non blocking receive.
//---------------------------------------------------------------------------//

C4_Req C4_RecvAsync( void *buf, int size, int source, int tag, int group )
{
    C4_Req r;

    if (source == C4_Any_Source)
	r.mid = irecv( tag, (char *) buf, size );
    else
	r.mid = irecvx( tag, (char *) buf, size, source, group, msginfo );

    r.set();
    return r;
}

//---------------------------------------------------------------------------//
// Optimized form which avoids spurious object creation and copy.
//---------------------------------------------------------------------------//

void C4_SendAsync( C4_Req& r, void *buf, int size, int dest, int tag, 
		   int group /*=0*/ )
{
// Not checking that r is not in use, which is of course a concern...
    r.mid = isend( tag, (char *) buf, size, dest, group );
    r.set();
}

//---------------------------------------------------------------------------//
// Optimized form which avoids spurious object creation and copy.
//---------------------------------------------------------------------------//

void C4_RecvAsync( C4_Req& r, void *buf, int size, int source, int tag, 
		   int group /*=0*/ )
{
// Not checking that r is not in use, which is of course a concern...
    if (source == C4_Any_Source)
	r.mid = irecv( tag, (char *) buf, size );
    else
	r.mid = irecvx( tag, (char *) buf, size, source, group, msginfo );

    r.set();
}

//---------------------------------------------------------------------------//
// Send a buffer of integers.
//---------------------------------------------------------------------------//

int C4_Send( int *buf, int nels, int dest, int group /*=0*/ )
{
    csend( C4_int_ptr_Tag, (char *) buf, nels*sizeof(int), dest, group );
    return C4_SUCCESS;
}

//---------------------------------------------------------------------------//
// Receive a buffer of integers.  nels is how many we can receive, the return
// value is how many we did receive.
//---------------------------------------------------------------------------//

int C4_Recv( int *buf, int nels, int source, int group /*=0*/ )
{
    int size = nels * sizeof(int);

    if (source == C4_Any_Source)
	crecv( C4_int_ptr_Tag, (char *) buf, size );
    else
	crecvx( C4_int_ptr_Tag, (char *) buf, size, source, group, msginfo );

    return infocount()/sizeof(int);
}

//---------------------------------------------------------------------------//
// Send a buffer of floats.
//---------------------------------------------------------------------------//

int C4_Send( float *buf, int nels, int dest, int group /*=0*/ )
{
    csend( C4_float_ptr_Tag, (char *) buf, nels*sizeof(float), dest, group );
    return C4_SUCCESS;
}

//---------------------------------------------------------------------------//
// Receive a buffer of floats.  nels is how many we can receive, the return
// value is how many we did receive.
//---------------------------------------------------------------------------//

int C4_Recv( float *buf, int nels, int source, int group /*=0*/ )
{
    int size = nels * sizeof(float);

    if (source == C4_Any_Source)
	crecv( C4_float_ptr_Tag, (char *) buf, size );
    else
	crecvx( C4_float_ptr_Tag, (char *) buf, size, source, group, msginfo );

    return infocount()/sizeof(float);
}

//---------------------------------------------------------------------------//
// Global reduction operations.
//
// Available operations:	description		types
//
// g<t>sum, g<t>prod		sum, product		i, s, d
// g<t>low, g<t>high		min, max		i, s, d
// giand, gior			bitwise and, or		(i only)
// gland, glor			logical and, or		(i only)
//
// where <t> is one of the listed types, corresponding to 
// i - long, s - float, d - double.
//
// All of the above functions take an argument list of the form:
// (x, n, work)
//
// Also there are operations: 
//	gcol, gcolx	concatenation
//	gopf		arbitrary commutative function
// (not described here)
//---------------------------------------------------------------------------//

static long   iwork;
static long   lwork;
static float  fwork;
static double dwork;

static DynArray<long>   ibuf(10);	// grrrrrrrrrr.
static DynArray<long>   lbuf(10);
static DynArray<float>  fbuf(10);
static DynArray<double> dbuf(10);

//---------------------------------------------------------------------------//
// Sum, scalar

void C4_gsum( int& x )
{
    long y = x;
    gisum( &y, 1, &iwork );
    x = (int) y;
}

void C4_gsum( long& x )
{
    gisum( &x, 1, &lwork );
}

void C4_gsum( float& x )
{
    gssum( &x, 1, &fwork );
}

void C4_gsum( double& x )
{
    gdsum( &x, 1, &dwork );
}

//---------------------------------------------------------------------------//
// Sum, array

void C4_gsum( int *px, int n )
{
    Assert( n >= 0 );
    lbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	ibuf[i] = px[i];

    gisum( &ibuf[0], n, &lbuf[0] );

    for( i=0; i < n; i++ )
	px[i] = (int) ibuf[i];
}

void C4_gsum( long *px, int n )
{
    Assert( n >= 0 );
    lbuf[n-1] = 0;		// auto expand the buffer.
    gisum( px, n, &lbuf[0] );
}

void C4_gsum( float *px, int n )
{
    Assert( n >= 0 );
    fbuf[n-1] = 0;		// auto expand the buffer.
    gssum( px, n, &fbuf[0] );
}

void C4_gsum( double *px, int n )
{
    Assert( n >= 0 );
    dbuf[n-1] = 0;		// auto expand the buffer.
    gdsum( px, n, &dbuf[0] );
}

//---------------------------------------------------------------------------//
// Min, scalar

void C4_gmin( int& x )
{
    long y = x;
    gilow( &y, 1, &iwork );
    x = (int) y;
}

void C4_gmin( long& x )
{
    gilow( &x, 1, &lwork );
}

void C4_gmin( float& x )
{
    gslow( &x, 1, &fwork );
}

void C4_gmin( double& x )
{
    gdlow( &x, 1, &dwork );
}

//---------------------------------------------------------------------------//
// Max, scalar

void C4_gmax( int& x )
{
    long y = x;
    gihigh( &y, 1, &iwork );
    x = (int) y;
}

void C4_gmax( long& x )
{
    gihigh( &x, 1, &lwork );
}

void C4_gmax( float& x )
{
    gshigh( &x, 1, &fwork );
}

void C4_gmax( double& x )
{
    gdhigh( &x, 1, &dwork );
}

//---------------------------------------------------------------------------//
//                              end of global_nx.cc
//---------------------------------------------------------------------------//
