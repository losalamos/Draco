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

const int C4_SUCCESS = 0;
#ifdef __PARAGON__
const int C4_Any_Tag = -1;
const int C4_Any_Source = -1;
#endif
#ifdef __C4_SCALAR__
const int C4_Any_Tag = -1;
const int C4_Any_Source = -1;
#endif
#ifdef C4_SHMEM
const int C4_Any_Tag = -1;
const int C4_Any_Source = -1;
#endif
#ifdef __MPI__
const int C4_Any_Tag = MPI_ANY_TAG;
const int C4_Any_Source = MPI_ANY_SOURCE;
#endif

const int C4_int_Tag = 432;
const int C4_float_Tag = 433;
const int C4_double_Tag = 434;
const int C4_int_ptr_Tag = 443;
const int C4_float_ptr_Tag = 444;
const int C4_double_ptr_Tag = 445;

enum { C4_Pass_Left, C4_Pass_Right };

//---------------------------------------------------------------------------//
// Prototypes
// Cruft

void C4_Init( int& argc, char **& argv );
void C4_Finalize();

// Informational 

int C4_node();
int C4_nodes();
int C4_group();

// Global sync

void C4_gsync();

// Send/receive

int C4_Send( void *buf, int size, int dest,   int tag, int group =0 );
int C4_Recv( void *buf, int size, int source, int tag, int group =0 );

C4_Req C4_SendAsync( void *buf, int size, int dest,   int tag, int group =0 );
C4_Req C4_RecvAsync( void *buf, int size, int source, int tag, int group =0 );

void C4_SendAsync( C4_Req& r, void *buf, int size, int dest,   int tag, int group =0 );
void C4_RecvAsync( C4_Req& r, void *buf, int size, int source, int tag, int group =0 );

// Super convenient special forms.

inline int C4_Send( int data, int dest, int group =0 )
{
    return C4_Send( &data, sizeof(int), dest, C4_int_Tag, group );
}

inline int C4_Recv( int& data, int source, int group =0 )
{
    return C4_Recv( &data, sizeof(int), source, C4_int_Tag, group );
}

inline int C4_Send( float data, int dest, int group =0 )
{
    return C4_Send( &data, sizeof(float), dest, C4_float_Tag, group );
}

inline int C4_Recv( float& data, int source, int group =0 )
{
    return C4_Recv( &data, sizeof(float), source, C4_float_Tag, group );
}

inline int C4_Send( double data, int dest, int group =0 )
{
    return C4_Send( &data, sizeof(double), dest, C4_double_Tag, group );
}

inline int C4_Recv( double& data, int source, int group =0 )
{
    return C4_Recv( &data, sizeof(double), source, C4_double_Tag, group );
}

int C4_Send( int *buf, int nels, int dest, int group =0 );
int C4_Recv( int *buf, int nels, int source, int group =0 );
int C4_Send( float *buf, int nels, int dest, int group =0 );
int C4_Recv( float *buf, int nels, int source, int group =0 );
int C4_Send( double *buf, int nels, int dest, int group =0 );
int C4_Recv( double *buf, int nels, int source, int group =0 );

// Global reductions
// Sum, scalar

void C4_gsum( int& x );
void C4_gsum( long& x );
void C4_gsum( float& x );
void C4_gsum( double& x );

// Sum, array

void C4_gsum( int *px, int n );
void C4_gsum( long *px, int n );
void C4_gsum( float *px, int n );
void C4_gsum( double *px, int n );

// Min, scalar

void C4_gmin( int& x );
void C4_gmin( long& x );
void C4_gmin( float& x );
void C4_gmin( double& x );

// Min, array

void C4_gmin( int *px, int n );
void C4_gmin( long *px, int n );
void C4_gmin( float *px, int n );
void C4_gmin( double *px, int n );

// Max, scalar

void C4_gmax( int& x );
void C4_gmax( long& x );
void C4_gmax( float& x );
void C4_gmax( double& x );

// Max, array

void C4_gmax( int *px, int n );
void C4_gmax( long *px, int n );
void C4_gmax( float *px, int n );
void C4_gmax( double *px, int n );

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
