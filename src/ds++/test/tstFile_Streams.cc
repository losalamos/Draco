//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstFile_Streams.cc
 * \author Rob Lowrie
 * \date   Sun Nov 21 19:36:12 2004
 * \brief  Tests File_Input and File_Output.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#include "../Assert.hh"
#include "../Release.hh"
#include "../File_Streams.hh"
#include "../Soft_Equivalence.hh"
#include "ds_test.hh"

using namespace std;
using rtt_dsxx::File_Input;
using rtt_dsxx::File_Output;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_fileio(const bool binary)
{
    string filename("file_streams.");

    if ( binary )
	filename += "binary";
    else
	filename += "ascii";

    int i = 5;
    string s = "  a string with spaces  ";
    double x = 5.6;
    bool bf = false; 
    bool bt = true;

    // write the data

    {
	File_Output f(filename, binary);
	f << i;

	// here's how you write strings:
	int ssize = s.size();
	f << ssize;
	for ( int k = 0; k < ssize; k++ )
	    f << s[k];

	f << x << bf << bt;
    }

    // read the data and make sure it's the same

    {
	int i_in;
	double x_in;
	string s_in;
	bool bf_in;
	bool bt_in;

	File_Input f(filename);
	f >> i_in;

	UNIT_TEST(i == i_in);

	// here's how you read strings:
	int ssize;
	f >> ssize;
	UNIT_TEST(ssize == int(s.size()));
	s_in.resize(ssize);
	for ( int k = 0; k < ssize; k++ )
	    f >> s_in[k];

	UNIT_TEST(s == s_in);

	f >> x_in >> bf_in >> bt_in;

	UNIT_TEST(soft_equiv(x, x_in));
	UNIT_TEST(bf == bf_in);
	UNIT_TEST(bt == bt_in);
    }

    if ( rtt_ds_test::passed )
    {
	ostringstream m;
	m << "test_fileio(";
	if ( binary )
	    m << "binary";
	else
	    m << "ascii";
	m << ") ok.";
	PASSMSG(m.str());
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    std::cout << argv[0] << ": version " 
		      << rtt_dsxx::release() 
		      << std::endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	test_fileio(false); // ascii
	test_fileio(true);  // binary
    }
    catch (std::exception &err)
    {
	std::cout << "ERROR: While testing tstFile_Streams, " 
		  << err.what()
		  << std::endl;
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing tstFile_Streams, " 
		  << "An unknown exception was thrown."
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
	      << std::endl;
    if (rtt_ds_test::passed) 
    {
        std::cout << "**** tstFile_Streams Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" 
	      << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstFile_Streams." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstFile_Streams.cc
//---------------------------------------------------------------------------//
