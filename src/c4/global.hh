//----------------------------------*-C++-*----------------------------------//
// global.hh
// Geoffrey Furnish
// Tue Dec 20 1994
//---------------------------------------------------------------------------//
// @> Improves on the NX API for global operations.
//---------------------------------------------------------------------------//

#ifndef __c4_global_hh__
#define __c4_global_hh__

#include "c4/config.hh"
#include "c4/C4_Req.hh"
#include "c4/tags.hh"

C4_NAMESPACE_BEG

void Init( int& argc, char **& argv );
void Finalize();

// Informational 

int node();
int nodes();
int group();

// Global sync

void gsync();

// Send/receive

int Send( void *buf, int size, int dest,   int tag, int group =0 );
int Recv( void *buf, int size, int source, int tag, int group =0 );

C4_Req SendAsync( void *buf, int size, int dest,   int tag, int group =0 );
C4_Req RecvAsync( void *buf, int size, int source, int tag, int group =0 );

void SendAsync( C4_Req& r, void *buf, int size, int dest,   int tag, int group =0 );
void RecvAsync( C4_Req& r, void *buf, int size, int source, int tag, int group =0 );

// Super convenient special forms.

inline int Send( int data, int dest, int group =0 )
{
    return Send( &data, sizeof(int), dest, C4_int_Tag, group );
}

inline int Recv( int& data, int source, int group =0 )
{
    return Recv( &data, sizeof(int), source, C4_int_Tag, group );
}

inline int Send( float data, int dest, int group =0 )
{
    return Send( &data, sizeof(float), dest, C4_float_Tag, group );
}

inline int Recv( float& data, int source, int group =0 )
{
    return Recv( &data, sizeof(float), source, C4_float_Tag, group );
}

inline int Send( double data, int dest, int group =0 )
{
    return Send( &data, sizeof(double), dest, C4_double_Tag, group );
}

inline int Recv( double& data, int source, int group =0 )
{
    return Recv( &data, sizeof(double), source, C4_double_Tag, group );
}
/*
int Send( int *buf, int nels, int dest, int group =0 );
int Recv( int *buf, int nels, int source, int group =0 );
int Send( float *buf, int nels, int dest, int group =0 );
int Recv( float *buf, int nels, int source, int group =0 );
int Send( double *buf, int nels, int dest, int group =0 );
int Recv( double *buf, int nels, int source, int group =0 );
*/

template<class T>
int Send( const T *buf, int nels, int dest, int group =0 );
template<class T>
int Recv( T *buf, int nels, int source, int group =0 );

// Global reductions
// Sum, scalar

void gsum( int& x );
void gsum( long& x );
void gsum( float& x );
void gsum( double& x );

// Sum, array

void gsum( int *px, int n );
void gsum( long *px, int n );
void gsum( float *px, int n );
void gsum( double *px, int n );

// Min, scalar

void gmin( int& x );
void gmin( long& x );
void gmin( float& x );
void gmin( double& x );

// Min, array

void gmin( int *px, int n );
void gmin( long *px, int n );
void gmin( float *px, int n );
void gmin( double *px, int n );

// Max, scalar

void gmax( int& x );
void gmax( long& x );
void gmax( float& x );
void gmax( double& x );

// Max, array

void gmax( int *px, int n );
void gmax( long *px, int n );
void gmax( float *px, int n );
void gmax( double *px, int n );

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
// Prototypes
// Cruft

inline void C4_Init( int& argc, char **& argv ) { C4::Init( argc, argv ); }
inline void C4_Finalize() { C4::Finalize(); }

// Informational 

inline int C4_node() { return C4::node(); }
inline int C4_nodes() { return C4::nodes(); }
inline int C4_group() { return C4::group(); }

// Global sync

inline void C4_gsync() { C4::gsync(); }

// Send/receive

inline int C4_Send( void *buf, int size, int dest,   int tag, int group =0 )
{ return C4::Send( buf, size, dest, tag, group ); }
inline int C4_Recv( void *buf, int size, int source, int tag, int group =0 )
{ return C4::Recv( buf, size, source, tag, group ); }

inline C4::C4_Req C4_SendAsync( void *buf, int size, int dest,
				int tag, int group =0 )
{ return C4::SendAsync( buf, size, dest, tag, group ); }
inline C4::C4_Req C4_RecvAsync( void *buf, int size, int source,
				int tag, int group =0 )
{ return C4::RecvAsync( buf, size, source, tag, group ); }

inline void C4_SendAsync( C4::C4_Req& r, void *buf, int size, int dest,
			  int tag, int group =0 )
{ C4::SendAsync( r, buf, size, dest, tag, group ); }
inline void C4_RecvAsync( C4::C4_Req& r, void *buf, int size, int source,
			  int tag, int group =0 )
{ C4::RecvAsync( r, buf, size, source, tag, group ); }

// Super convenient special forms.

inline int C4_Send( int data, int dest, int group =0 )
{
    return C4::Send( &data, sizeof(int), dest, C4_int_Tag, group );
}

inline int C4_Recv( int& data, int source, int group =0 )
{
    return C4::Recv( &data, sizeof(int), source, C4_int_Tag, group );
}

inline int C4_Send( float data, int dest, int group =0 )
{
    return C4::Send( &data, sizeof(float), dest, C4_float_Tag, group );
}

inline int C4_Recv( float& data, int source, int group =0 )
{
    return C4::Recv( &data, sizeof(float), source, C4_float_Tag, group );
}

inline int C4_Send( double data, int dest, int group =0 )
{
    return C4::Send( &data, sizeof(double), dest, C4_double_Tag, group );
}

inline int C4_Recv( double& data, int source, int group =0 )
{
    return C4::Recv( &data, sizeof(double), source, C4_double_Tag, group );
}

inline int C4_Send( int *buf, int nels, int dest, int group =0 )
{ return C4::Send( buf, nels, dest, group ); }
inline int C4_Recv( int *buf, int nels, int source, int group =0 )
{ return C4::Recv( buf, nels, source, group ); }
inline int C4_Send( float *buf, int nels, int dest, int group =0 )
{ return C4::Send( buf, nels, dest, group ); }
inline int C4_Recv( float *buf, int nels, int source, int group =0 )
{ return C4::Recv( buf, nels, source, group ); }
inline int C4_Send( double *buf, int nels, int dest, int group =0 )
{ return C4::Send( buf, nels, dest, group ); }
inline int C4_Recv( double *buf, int nels, int source, int group =0 )
{ return C4::Recv( buf, nels, source, group ); }

// Global reductions
// Sum, scalar

inline void C4_gsum( int& x ) { C4::gsum( x ); }
inline void C4_gsum( long& x ) { C4::gsum( x ); }
inline void C4_gsum( float& x ) { C4::gsum( x ); }
inline void C4_gsum( double& x ) { C4::gsum( x ); }

// Sum, array

inline void C4_gsum( int *px, int n ) { C4::gsum( px, n ); }
inline void C4_gsum( long *px, int n ) { C4::gsum( px, n ); }
inline void C4_gsum( float *px, int n ) { C4::gsum( px, n ); }
inline void C4_gsum( double *px, int n ) { C4::gsum( px, n ); }

// Min, scalar

inline void C4_gmin( int& x ) { C4::gmin( x ); }
inline void C4_gmin( long& x ) { C4::gmin( x ); }
inline void C4_gmin( float& x ) { C4::gmin( x ); }
inline void C4_gmin( double& x ) { C4::gmin( x ); }

// Min, array

inline void C4_gmin( int *px, int n ) { C4::gmin( px, n ); }
inline void C4_gmin( long *px, int n ) { C4::gmin( px, n ); }
inline void C4_gmin( float *px, int n ) { C4::gmin( px, n ); }
inline void C4_gmin( double *px, int n ) { C4::gmin( px, n ); }

// Max, scalar

inline void C4_gmax( int& x ) { C4::gmax( x ); }
inline void C4_gmax( long& x ) { C4::gmax( x ); }
inline void C4_gmax( float& x ) { C4::gmax( x ); }
inline void C4_gmax( double& x ) { C4::gmax( x ); }

// Max, array

inline void C4_gmax( int *px, int n ) { C4::gmax( px, n ); }
inline void C4_gmax( long *px, int n ) { C4::gmax( px, n ); }
inline void C4_gmax( float *px, int n ) { C4::gmax( px, n ); }
inline void C4_gmax( double *px, int n ) { C4::gmax( px, n ); }

//---------------------------------------------------------------------------//
// Now include anything which might help a particular hardware abstraction
// layer.  For example, C4_gsync might be inlined...

#ifdef C4_SHMEM
#include "c4/global_shmem.hh"
#endif

#endif                          // __c4_global_hh__

//---------------------------------------------------------------------------//
//                              end of c4/global.hh
//---------------------------------------------------------------------------//
