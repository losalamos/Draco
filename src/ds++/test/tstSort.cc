//----------------------------------*-C++-*----------------------------------//
// tstSort.cc
// Geoffrey Furnish
// 28 September 1995
//---------------------------------------------------------------------------//
// @> Program to test the sorting functions.
//---------------------------------------------------------------------------//

#include "../Sort.hh"

#include <iostream>

#include <cstdlib>
#include <cstdarg>

using namespace std;

void test_HeapSort()
{
    float *v = new float[ 100 ];
    int i;

    for( i=0; i < 10; i++ )
	v[i] = 10. * rand() / RAND_MAX;

    HeapSort( 10, v );

    for( i=0; i < 10; i++ )
	cout << i << ' ' << v[i] << endl;
}

void test_index()
{
    cout << "Testing index\n";

    int *v = new int[ 10 ];
    int i;

    for( i=0; i < 10; i++ )
	v[i] = static_cast<int>(10. * rand() / RAND_MAX); 

    cout << "Starting data is:\n";
    for( i=0; i < 10; i++ )
	cout << i << ' ' << v[i] << endl;

    int *ndx = new int[ 10 ];

    index( 10, v, ndx );

    cout << "Indexed data is:\n";
    for( i=0; i < 10; i++ )
	cout << i << ' ' << v[ ndx[i] ] << endl;

    cout << "Done testing index.\n";
}    

int main( int argc, char *argv[] )
{
    cout << "Sort test:\n";

//     test_HeapSort();
    test_index();

    cout << "Done.\n";

    cout << "tstSort Test: passed" << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstSort.cc
//---------------------------------------------------------------------------//
