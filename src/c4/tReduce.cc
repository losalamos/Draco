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

//---------------------------------------------------------------------------//

template<class T>
void tReduce_sum( T dummy )
{
    T local_sum = 0, global_sum = 0;

    for (int i = 0; i < nodes; i++)
	local_sum += i;

    global_sum = node;
    C4::gsum(global_sum);

    if (node == 0) {
	if (local_sum == global_sum)
	    cout << "gsum for type "
		 << typeid(T).name()
		 << " correct: " << global_sum << endl;
	else
	    cout << "gsum for type "
		 << typeid(T).name()
		 << " bogus: " << global_sum << endl;
    }
}

//---------------------------------------------------------------------------//

template<class T>
void tReduce_min( T dummy )
{
    int local_min = -nodes+1, global_min;

    global_min = -node;
    C4::gmin(global_min);

    if (node == 0) {
	if (local_min == global_min)
	    cout << "gmin for type "
		 << typeid(T).name()
		 << " correct: " << global_min << endl;
	else
	    cout << "gmin for type "
		 << typeid(T).name()
		 << " bogus: " << global_min << endl;
    }
}

//---------------------------------------------------------------------------//

template<class T>
void tReduce_max( T dummy )
{
    int local_max = nodes-1, global_max;

    global_max = node;
    C4::gmax(global_max);

    if (node == 0) {
	if (local_max == global_max)
	    cout << "gmax for type "
		 << typeid(T).name()
		 << " correct: " << global_max << endl;
	else
	    cout << "gmax for type "
		 << typeid(T).name()
		 << " bogus: " << global_max << endl;
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

    tReduce_sum( int() );
    tReduce_sum( float() );
    tReduce_sum( double() );

    tReduce_min( int() );
    tReduce_min( float() );
    tReduce_min( double() );

    tReduce_max( int() );
    tReduce_max( float() );
    tReduce_max( double() );

    C4::Finalize();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tReduce.cc
//---------------------------------------------------------------------------//
