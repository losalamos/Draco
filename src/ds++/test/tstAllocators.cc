//----------------------------------*-C++-*----------------------------------//
// tstAllocators.cc
// Geoffrey M. Furnish
// Mon Mar  2 11:09:12 1998
//---------------------------------------------------------------------------//
// @> Exercise the DS++ Allocators.hh component.
//---------------------------------------------------------------------------//

#include "../Allocators.hh"

#include <iostream>
#include <string>

using namespace rtt_dsxx;
using namespace std;

//---------------------------------------------------------------------------//
// Register a test as passed.
//---------------------------------------------------------------------------//

template<class T>
void pass( const char *name, T dummy )
{
    cout << name << '<' << typeid(T).name() << ">test: passed" << endl;
}

//---------------------------------------------------------------------------//
// Register a test as failed.
//---------------------------------------------------------------------------//

template<class T>
void fail( const char *name, T dummy )
{
    cout << name << '<' << typeid(T).name() << ">test: failed" << endl;
}

//---------------------------------------------------------------------------//
// Test basic operation of Simple_Allocator.
//---------------------------------------------------------------------------//

template<class T>
void tS1( T dummy )
{
    try {
        T *p = Simple_Allocator<T>::fetch( 5 );
        Simple_Allocator<T>::release( p, 5 );
	pass( "tS1", dummy );
    }
    catch(...)
    {
	fail( "tS1", dummy );
    }
}

//---------------------------------------------------------------------------//
// Test basic operation of Guarded_Allocator when nothing bogus is going on.
//---------------------------------------------------------------------------//

template<class T>
void tG1( T dummy )
{
    try {
        T *p = Guarded_Allocator<T>::fetch( 5 );
        Guarded_Allocator<T>::release( p, 5 );
	pass( "tG1", dummy );
    }
    catch(...)
    {
	fail( "tG1", dummy );
    }
}

//---------------------------------------------------------------------------//
// Check that Guarded_Allocator can detect subversive behavior.
//---------------------------------------------------------------------------//

template<class T>
void tG2( T dummy )
{
    try {
    // First we have to fetch some memory.
        T *p = Guarded_Allocator<T>::fetch( 5 );

    // Now initialize it.
	std::uninitialized_fill_n( p, 5, dummy );

    // Now currupt the bottom end of the memory :-).
	p[-1] = dummy;

    // Now release the memory.
        Guarded_Allocator<T>::release( p, 5 );
#if DBC & 2
	fail( "tG2", dummy );
#else
	pass( "tG2", dummy );
#endif
    }
    catch(assertion& x)
    {
#if DBC & 2
	pass( "tG2", dummy );
#else
	fail( "tG2", dummy );
#endif
    }
    catch(...)
    {
	fail( "tG2", dummy );
    }
}

//---------------------------------------------------------------------------//
// Check that Guarded_Allocator can detect subversive behavior.
//---------------------------------------------------------------------------//

template<class T>
void tG3( T dummy )
{
    try {
    // First we have to fetch some memory.
        T *p = Guarded_Allocator<T>::fetch( 5 );

    // Now initialize it.
	std::uninitialized_fill_n( p, 5, dummy );

    // Now currupt the top end of the memory :-).
	p[5] = dummy;

    // Now release the memory.
        Guarded_Allocator<T>::release( p, 5 );
#if DBC & 2
	fail( "tG3", dummy );
#else
	pass( "tG3", dummy );
#endif
    }
    catch(assertion& x)
    {
#if DBC & 2
	pass( "tG3", dummy );
#else
	fail( "tG3", dummy );
#endif
    }
    catch(...)
    {
	fail( "tG3", dummy );
    }
}

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main( int argc, char *argv[] )
{

    for (int arg=1; arg < argc; arg++)
	{
	    if (std::string(argv[arg]) == "--version")
		{
		    version(argv[0]);
		    return 0;
		}
	}

    cout << "Starting tstAllocators.\n";

    tS1( 7 );

    tG1( 7 );
    tG2( 7 );
    tG3( 7 );

    cout << "Ending tstAllocators.\n";
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstAllocators.cc
//---------------------------------------------------------------------------//
