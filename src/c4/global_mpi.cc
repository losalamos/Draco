//----------------------------------*-C++-*----------------------------------//
// global_mpi.cc
// Maurice LeBrun
// Wed Feb  1 16:50:59 1995
//---------------------------------------------------------------------------//
// @> Global C4 functions for MPI
//---------------------------------------------------------------------------//

#include "ds++/DynArray.hh"
#include "ds++/Assert.hh"

//---------------------------------------------------------------------------//
// Miscellaneous

C4_NAMESPACE_BEG

void Init( int& argc, char **& argv )
{
    MPI_Init( &argc, &argv );
}

void Finalize()
{
    MPI_Finalize();
}

int node()
{
    int node;
    MPI_Comm_rank( MPI_COMM_WORLD, &node );
    return node;
}

int nodes()
{
    int nodes;
    MPI_Comm_size( MPI_COMM_WORLD, &nodes );
    return nodes;
}

int group()
{
    int group = 0;
    return group;
}

void gsync()
{
    MPI_Barrier( MPI_COMM_WORLD );
}

//---------------------------------------------------------------------------//
// Global reduction operations.
//
// The call most like that in NX is:
//
//	MPI_Allreduce(void* sendbuf, void* recvbuf, int count,
//		      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
// 
// which returns the result to all processes in the group.  "op" determines
// the type of reduction performed.
//
// Available reduction operators, description, and allowed types: 
//
// MPI_MAX, MPI_MIN		min, max			I, F
// MPI_SUM, MPI_PROD		sum, product			I, F
// MPI_BAND, MPI_BOR, MPI_BXOR	bitwise and, or, xor		I, B
// MPI_LAND, MPI_LOR, MPI_LXOR	logical and, or, xor		I
// MPI_MAXLOC, MPI_MINLOC	min, max value and location
//
// where types are:
//
// I:	MPI_INT, MPI_LONG, MPI_SHORT, MPI_UNSIGNED_SHORT, 
//	MPI_UNSIGNED, MPI_UNSIGNED_LONG
// F:	MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE
// B:	MPI_BYTE
//---------------------------------------------------------------------------//

static DynArray<int>    ibuf(10);
static DynArray<long>   lbuf(10);
static DynArray<float>  fbuf(10);
static DynArray<double> dbuf(10);

//---------------------------------------------------------------------------//
// Sum, scalar

void gsum( int& x )
{
    int y = x;
    MPI_Allreduce( &y, &x, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
}

void gsum( long& x )
{
    long y = x;
    MPI_Allreduce( &y, &x, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
}

void gsum( float& x )
{
    float y = x;
    MPI_Allreduce( &y, &x, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
}

void gsum( double& x )
{
    double y = x;
    MPI_Allreduce( &y, &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
}

//---------------------------------------------------------------------------//
// Sum, array

void gsum( int *px, int n )
{
    Assert( n >= 0 );

    ibuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	ibuf[i] = px[i];

    MPI_Allreduce( &ibuf[0], px, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
}

void gsum( long *px, int n )
{
    Assert( n >= 0 );

    lbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	lbuf[i] = px[i];

    MPI_Allreduce( &lbuf[0], px, n, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
}

void gsum( float *px, int n )
{
    Assert( n >= 0 );

    fbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	fbuf[i] = px[i];

    MPI_Allreduce( &fbuf[0], px, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );
}

void gsum( double *px, int n )
{
    Assert( n >= 0 );

    dbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	dbuf[i] = px[i];

    MPI_Allreduce( &dbuf[0], px, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
}

//---------------------------------------------------------------------------//
// Min, scalar

void gmin( int& x )
{
    int y = x;
    MPI_Allreduce( &y, &x, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
}

void gmin( long& x )
{
    long y = x;
    MPI_Allreduce( &y, &x, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD );
}

void gmin( float& x )
{
    float y = x;
    MPI_Allreduce( &y, &x, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD );
}

void gmin( double& x )
{
    double y = x;
    MPI_Allreduce( &y, &x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
}

//---------------------------------------------------------------------------//
// Min, array

void gmin( int *px, int n )
{
    Assert( n >= 0 );

    ibuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	ibuf[i] = px[i];

    MPI_Allreduce( &ibuf[0], px, n, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
}

void gmin( long *px, int n )
{
    Assert( n >= 0 );

    lbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	lbuf[i] = px[i];

    MPI_Allreduce( &lbuf[0], px, n, MPI_LONG, MPI_MIN, MPI_COMM_WORLD );
}

void gmin( float *px, int n )
{
    Assert( n >= 0 );

    fbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	fbuf[i] = px[i];

    MPI_Allreduce( &fbuf[0], px, n, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD );
}

void gmin( double *px, int n )
{
    Assert( n >= 0 );

    dbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	dbuf[i] = px[i];

    MPI_Allreduce( &dbuf[0], px, n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
}

//---------------------------------------------------------------------------//
// Max, scalar

void gmax( int& x )
{
    int y = x;
    MPI_Allreduce( &y, &x, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
}

void gmax( long& x )
{
    long y = x;
    MPI_Allreduce( &y, &x, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD );
}

void gmax( float& x )
{
    float y = x;
    MPI_Allreduce( &y, &x, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );
}

void gmax( double& x )
{
    double y = x;
    MPI_Allreduce( &y, &x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
}

//---------------------------------------------------------------------------//
// Max, array

void gmax( int *px, int n )
{
    Assert( n >= 0 );

    ibuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	ibuf[i] = px[i];

    MPI_Allreduce( &ibuf[0], px, n, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
}

void gmax( long *px, int n )
{
    Assert( n >= 0 );

    lbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	lbuf[i] = px[i];

    MPI_Allreduce( &lbuf[0], px, n, MPI_LONG, MPI_MAX, MPI_COMM_WORLD );
}

void gmax( float *px, int n )
{
    Assert( n >= 0 );

    fbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	fbuf[i] = px[i];

    MPI_Allreduce( &fbuf[0], px, n, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );
}

void gmax( double *px, int n )
{
    Assert( n >= 0 );

    dbuf[n-1] = 0;		// auto expand the buffer.
    for( int i=0; i < n; i++ )
	dbuf[i] = px[i];

    MPI_Allreduce( &dbuf[0], px, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of global_mpi.cc
//---------------------------------------------------------------------------//
