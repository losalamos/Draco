//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/test/tstEnsight_Stream.cc
 * \author Rob Lowrie
 * \date   Fri Nov 12 22:52:46 2004
 * \brief  Test for Ensight_Stream.
 * \note   Copyright 2004-2006 The Regents of the University of California.
 *         Copyright 2006-2010 LANS, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"
#include "viz_test.hh"

#include "../Ensight_Stream.hh"

using namespace std;
using rtt_viz::Ensight_Stream;

//---------------------------------------------------------------------------//
// Utility functions
//---------------------------------------------------------------------------//

// Reads binary value from stream.
template <class T>
void binary_read(ifstream &stream,
		 T &v)
{
    char *vc = new char[sizeof(T)];
    stream.read(vc, sizeof(T));
    
    rtt_dsxx::Unpacker p;
    p.set_buffer(sizeof(T), vc);
    p.unpack(v);

    delete[] vc;
}

// Various overloaded read functions.

void readit(ifstream &stream,
	    const bool binary,
	    double &d)
{
    if ( binary )
    {
	float x;
	binary_read(stream, x);
	d = x;
    }
    else
	stream >> d;
}

void readit(ifstream &stream,
	    const bool binary,
	    int &d)
{
    if ( binary )
	binary_read(stream, d);
    else
	stream >> d;
}

void readit(ifstream &stream,
	    const bool binary,
	    string &s)
{
    if ( binary )
    {
	s.resize(80);
        for ( int i = 0; i < 80; ++ i )
	    stream.read(&s[i], 1);
    }
    else
	stream >> s;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_simple(const bool binary)
{
    // Dump a few values into the stream

    const int i = 20323;
    const string s("dog");
    const double d = 112.3;
    const string file("ensight_stream.out");
    
    {
	Ensight_Stream f(file, binary);
	
	f << i << rtt_viz::endl;
	f << d << rtt_viz::endl;
	f << s << rtt_viz::endl;
    }

    // Read the file back in and check the values.
    
    std::ios::openmode mode = std::ios::in;

    if ( binary )
    {
	cout << "Testing binary mode." << endl;
        mode = mode | std::ios::binary;
    }
    else
	cout << "Testing ascii mode." << endl;

    ifstream in(file.c_str(), mode);
    
    int i_in;
    readit(in, binary, i_in);
    UNIT_TEST(i == i_in);
    
    double d_in;
    readit(in, binary, d_in);
    UNIT_TEST(rtt_dsxx::soft_equiv(d, d_in, 0.01)); // floats are inaccurate

    string s_in;
    readit(in, binary, s_in);
    for ( size_t k = 0; k < s.size(); ++k )
	UNIT_TEST(s[k] == s_in[k]);
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    std::cout << argv[0] << ": version " << rtt_dsxx::release() 
              << std::endl;
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	    return 0;
    
    try
    {
	// >>> UNIT TESTS

	test_simple(true);  // test binary
	test_simple(false); // test ascii
	test_simple(true);  // test binary again
    }
    catch (std::exception &err)
    {
	std::cout << "ERROR: While testing tstEnsight_Stream, " 
		  << err.what()
		  << std::endl;
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing tstEnsight_Stream, " 
		  << "An unknown exception was thrown."
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
	      << std::endl;
    if (rtt_viz_test::passed) 
    {
        std::cout << "**** tstEnsight_Stream Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" 
	      << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstEnsight_Stream." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstEnsight_Stream.cc
//---------------------------------------------------------------------------//
