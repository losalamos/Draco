//----------------------------------*-C++-*----------------------------------//
// tReduce.cc
// Maurice LeBrun
// Mon Feb  6 15:16:26 1995
//---------------------------------------------------------------------------//
// @> Test program for reduction operations.
//---------------------------------------------------------------------------//

#include "c4/global.hh"

#include <iostream.h>

int node, nodes;

template<class T>
bool equal_buffers( T *s1, const T *e1, T *s2 )
{
    while( s1 < e1 )
	if (*s1++ != *s2++) return false;

    return true;
}

//---------------------------------------------------------------------------//

template<class T>
void scalar_sum( T dummy )
{
    T local_sum = 0, global_sum = 0;

    for (int i = 0; i < nodes; i++)
	local_sum += i;

    global_sum = node;
    C4::gsum(global_sum);

    if (node == 0) {
	if (local_sum == global_sum)
	    cout << "scalar sum for type "
		 << typeid(T).name()
		 << " correct: " << global_sum << endl;
	else
	    cout << "scalar sum for type "
		 << typeid(T).name()
		 << " bogus: " << global_sum << endl;
    }
}

//---------------------------------------------------------------------------//

template<class T>
void scalar_min( T dummy )
{
    int local_min = -nodes+1, global_min;

    global_min = -node;
    C4::gmin(global_min);

    if (node == 0) {
	if (local_min == global_min)
	    cout << "scalar min for type "
		 << typeid(T).name()
		 << " correct: " << global_min << endl;
	else
	    cout << "scalar min for type "
		 << typeid(T).name()
		 << " bogus: " << global_min << endl;
    }
}

//---------------------------------------------------------------------------//

template<class T>
void scalar_max( T dummy )
{
    int local_max = nodes-1, global_max;

    global_max = node;
    C4::gmax(global_max);

    if (node == 0) {
	if (local_max == global_max)
	    cout << "scalar max for type "
		 << typeid(T).name()
		 << " correct: " << global_max << endl;
	else
	    cout << "scalar max for type "
		 << typeid(T).name()
		 << " bogus: " << global_max << endl;
    }
}

//---------------------------------------------------------------------------//

template<class T>
void array_sum( T dummy )
{
    T *global_sum = new T[ 20 ];
    T *local_sum = new T[ 20 ];

    for( int i=0; i < 20; i++ ) {
	local_sum[i] = 0;
	global_sum[i] = 100*node + i;
	for( int j=0; j < nodes; j++ )
	    local_sum[i] +=  100*j + i;
    }

    C4::gsum( global_sum, 20 );

    if (node == 0) {
	if ( equal_buffers( local_sum, local_sum+20, global_sum ) )
	    cout << "array sum for type "
		 << typeid(T).name()
		 << ": correct." << endl;
	else
	    cout << "array sum for type "
		 << typeid(T).name()
		 << ": bogus." << endl;
    }
}

//---------------------------------------------------------------------------//

template<class T>
void array_min( T dummy )
{
    T *global_min = new T[ 20 ];
    T *local_min = new T[ 20 ];

    for( int i=0; i < 20; i++ ) {
	global_min[i] = -node + i;
	local_min[i] = -nodes+1 + i;
    }

    C4::gmin( global_min, 20 );

    if (node == 0) {
	if ( equal_buffers( local_min, local_min+20, global_min ) )
	    cout << "array min for type "
		 << typeid(T).name()
		 << ": correct." << endl;
	else
	    cout << "array min for type "
		 << typeid(T).name()
		 << ": bogus." << endl;
    }
}

//---------------------------------------------------------------------------//

template<class T>
void array_max( T dummy )
{
    T *global_max = new T[ 20 ];
    T *local_max = new T[ 20 ];

    for( int i=0; i < 20; i++ ) {
	global_max[i] = 100*i + node;
	local_max[i] = 100*i + nodes - 1;
    }

    C4::gmax( global_max, 20 );

    if (node == 0) {
	if ( equal_buffers( local_max, local_max+20, global_max ) )
	    cout << "array max for type "
		 << typeid(T).name()
		 << ": correct." << endl;
	else
	    cout << "array max for type "
		 << typeid(T).name()
		 << ": bogus." << endl;
    }
}

//---------------------------------------------------------------------------//
// Verify that global reduction operations work.
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    node = C4::node();
    nodes = C4::nodes();

    scalar_sum( int() );
    scalar_sum( float() );
    scalar_sum( double() );

    scalar_min( int() );
    scalar_min( float() );
    scalar_min( double() );

    scalar_max( int() );
    scalar_max( float() );
    scalar_max( double() );

    array_sum( int() );
    array_sum( float() );
    array_sum( double() );

    array_min( int() );
    array_min( float() );
    array_min( double() );

    array_max( int() );
    array_max( float() );
    array_max( double() );

    C4::Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tReduce.cc
//---------------------------------------------------------------------------//
