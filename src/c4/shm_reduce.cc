//----------------------------------*-C++-*----------------------------------//
// shm_reduce.cc
// Geoffrey M. Furnish
// Wed Feb 25 08:41:12 1998
//---------------------------------------------------------------------------//
// @> Implement the C4 global reduction API over SHMEM.
//---------------------------------------------------------------------------//

#include "c4/shmem.hh"
#include "c4/global.hh"

#include <mpp/shmem.h>

#include <algorithm>

C4_NAMESPACE_BEG

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

// "sr" == "scalar reduce"

static int sri_pWrk[ _SHMEM_REDUCE_MIN_WRKDATA_SIZE ];
static float srf_pWrk[ _SHMEM_REDUCE_MIN_WRKDATA_SIZE ];
static double srd_pWrk[ _SHMEM_REDUCE_MIN_WRKDATA_SIZE ];
static long sr_pSync[ _SHMEM_REDUCE_SYNC_SIZE ];

void C4_shm_init_scalar_work_arrays()
{
    for( int i=0; i < _SHMEM_REDUCE_SYNC_SIZE; i++ ) {
	sr_pSync[i] = _SHMEM_SYNC_VALUE;
    }

    gsync();
}

const int ar_data_size = 256;
#if ar_data_suze/2 + 1 > _SHMEM_REDUCE_MIN_WRKDATA_SIZE
const int ar_pWrk_size = ar_data_size/2 +1
#else
const int ar_pWrk_size = _SHMEM_REDUCE_MIN_WRKDATA_SIZE;
#endif

template<class T> class shm_ar_helper
{
  protected:
    static T ar_data[ ar_data_size ];
    static T ar_pWrk[ ar_pWrk_size ];

  public:
    static void copy_data_to_symmetric_work_arrays( T *px, int n )
    {
	Assert( n <= ar_data_size );
	for( int i=0; i < n; i++ )
	    ar_data[i] = px[i];
    }

    static void get_data_from_symmetric_work_arrays( T *px, int n )
    {
	Assert( n <= ar_data_size );
	for( int i=0; i < n; i++ )
	    px[i] = ar_data[i];
    }
};

template<class T>
T shm_ar_helper<T>::ar_data[ar_data_size];
template<class T>
T shm_ar_helper<T>::ar_pWrk[ar_pWrk_size];

template<class T> class shm_traits {};

template<> class shm_traits<int> : public shm_ar_helper<int>
{
  public:
    static void sum_to_all( int& xo, int& xi )
    {
	shmem_int_sum_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
			      sri_pWrk, sr_pSync );
	gsync();
    }
    static void min_to_all( int& xo, int& xi )
    {
	shmem_int_min_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
			      sri_pWrk, sr_pSync );
	gsync();
    }
    static void max_to_all( int& xo, int& xi )
    {
	shmem_int_max_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
			      sri_pWrk, sr_pSync );
	gsync();
    }
    static void ar_sum_to_all( int n )
    {
	shmem_int_sum_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
			      ar_pWrk, sr_pSync );
	gsync();
    }
    static void ar_min_to_all( int n )
    {
	shmem_int_min_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
			      ar_pWrk, sr_pSync );
	gsync();
    }
    static void ar_max_to_all( int n )
    {
	shmem_int_max_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
			      ar_pWrk, sr_pSync );
	gsync();
    }
};

template<> class shm_traits<float> : public shm_ar_helper<float>
{
  public:
    static void sum_to_all( float& xo, float& xi )
    {
	shmem_float_sum_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
				srf_pWrk, sr_pSync );
	gsync();
    }
    static void min_to_all( float& xo, float& xi )
    {
	shmem_float_min_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
				srf_pWrk, sr_pSync );
	gsync();
    }
    static void max_to_all( float& xo, float& xi )
    {
	shmem_float_max_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
				srf_pWrk, sr_pSync );
	gsync();
    }
    static void ar_sum_to_all( int n )
    {
	shmem_float_sum_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
				ar_pWrk, sr_pSync );
	gsync();
    }
    static void ar_min_to_all( int n )
    {
	shmem_float_min_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
				ar_pWrk, sr_pSync );
	gsync();
    }
    static void ar_max_to_all( int n )
    {
	shmem_float_max_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
				ar_pWrk, sr_pSync );
	gsync();
    }
};

template<> class shm_traits<double> : public shm_ar_helper<double>
{
  public:
    static void sum_to_all( double& xo, double& xi )
    {
	shmem_double_sum_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
				 srd_pWrk, sr_pSync );
	gsync();
    }
    static void min_to_all( double& xo, double& xi )
    {
	shmem_double_min_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
				 srd_pWrk, sr_pSync );
	gsync();
    }
    static void max_to_all( double& xo, double& xi )
    {
	shmem_double_max_to_all( &xo, &xi, 1, 0, 0, C4_shm_nodes,
				 srd_pWrk, sr_pSync );
	gsync();
    }
    static void ar_sum_to_all( int n )
    {
	shmem_double_sum_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
				 ar_pWrk, sr_pSync );
	gsync();
    }
    static void ar_min_to_all( int n )
    {
	shmem_double_min_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
				 ar_pWrk, sr_pSync );
	gsync();
    }
    static void ar_max_to_all( int n )
    {
	shmem_double_max_to_all( ar_data, ar_data, n, 0, 0, C4_shm_nodes,
				 ar_pWrk, sr_pSync );
	gsync();
    }
};

//---------------------------------------------------------------------------//
// Scalar reduction API.

template<class T>
void gsum( T& x )
{
    static T s;
    s = x;
    shm_traits<T>::sum_to_all( s, s );
    x = s;
}

template<class T>
void gmin( T& x )
{
    static T s;
    s = x;
    shm_traits<T>::min_to_all( s, s );
    x = s;
}

template<class T>
void gmax( T& x )
{
    static T s;
    s = x;
    shm_traits<T>::max_to_all( s, s );
    x = s;
}

template void gsum( int& x );
template void gsum( float& x );
template void gsum( double& x );

template void gmin( int& x );
template void gmin( float& x );
template void gmin( double& x );

template void gmax( int& x );
template void gmax( float& x );
template void gmax( double& x );

//---------------------------------------------------------------------------//
// Array reduction API.
//
// Because the sucky SGI SHMEM implementation doesn't have shmalloc, we have
// to "packetize" these reduction ops.  Grrrr.
//---------------------------------------------------------------------------//

template<class T> 
void gsum( T *px, int n )
{
    for( int i=0; i < n; i += ar_data_size )
    {
	int ilim = std::min( ar_data_size, n - i );
	shm_traits<T>::copy_data_to_symmetric_work_arrays( px+i, ilim );
	shm_traits<T>::ar_sum_to_all( ilim );
	shm_traits<T>::get_data_from_symmetric_work_arrays( px+i, ilim );
    }
}

template<class T> 
void gmin( T *px, int n )
{
    for( int i=0; i < n; i += ar_data_size )
    {
	int ilim = std::min( ar_data_size, n - i );
	shm_traits<T>::copy_data_to_symmetric_work_arrays( px+i, ilim );
	shm_traits<T>::ar_min_to_all( ilim );
	shm_traits<T>::get_data_from_symmetric_work_arrays( px+i, ilim );
    }
}

template<class T> 
void gmax( T *px, int n )
{
    for( int i=0; i < n; i += ar_data_size )
    {
	int ilim = std::min( ar_data_size, n - i );
	shm_traits<T>::copy_data_to_symmetric_work_arrays( px+i, ilim );
	shm_traits<T>::ar_max_to_all( ilim );
	shm_traits<T>::get_data_from_symmetric_work_arrays( px+i, ilim );
    }
}

template void gsum( int *, int );
template void gsum( float *, int );
template void gsum( double *, int );

template void gmin( int *, int );
template void gmin( float *, int );
template void gmin( double *, int );

template void gmax( int *, int );
template void gmax( float *, int );
template void gmax( double *, int );

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of shm_reduce.cc
//---------------------------------------------------------------------------//
