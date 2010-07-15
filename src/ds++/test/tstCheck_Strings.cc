//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstCheck_Strings.cc
 * \author John M. McGhee
 * \date   Sun Jan 30 14:57:09 2000
 * \brief  Test code for the Check_Strings utility functions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Check_Strings.hh"
#include "../Release.hh"
#include "../Assert.hh"
#include <iostream>
#include <string>
#include <vector>

//---------------------------------------------------------------------------//

static void Check_Strings_test()
{
    // Define a vector of strings for testing
    std::string n[] = {"this", "is", "a#", "test", "xxx!", "space check",
		       "123", "x", "test", "dog", "is", "cat", "", "abc" };
    const int nn = sizeof(n)/sizeof(std::string);
    std::vector<std::string> names(&n[0], &n[nn]);
    typedef std::vector<std::string>::iterator VS_iter;
    bool passed;

    // Print a header
    std::cout << std::endl << "*** String Utilities Test Program ***" 
	      << std::endl << std::endl;

    // List the test string

    std::cout << "The " << names.size() << 
	" strings to be tested are: " << std::endl;
    for (size_t i=0; i<names.size(); ++i)
	std::cout  << "\"" << names[i]  << "\"" << std::endl;
	
    std::cout << std::endl;

    //---------------------------------------------------------------------------//

    // Test for illegal characters.

    std::cout << "Illegal character utility test:" << std::endl;
    std::string bad_chars = "()[]* !^#$/";
    std::vector<VS_iter> result =
	rtt_dsxx::check_string_chars(names.begin(), names.end(), bad_chars);
    if (result.size() == 0) 
	std::cout << "All characters OK!" << std::endl;
    else
    {
	std::cout << "*** Error in string definition -" << std::endl;
	for (size_t i=0; i<result.size(); i++)
	    std::cout << "Found disallowed character(s) in string: \"" 
		      << *result[i] << "\"" << std::endl;
	std::cout << "The following characters are forbidden:" << 
	    std::endl << " \"" << bad_chars << "\"," << 
	    " as well as any white-space characters." << std::endl;
    }

    passed = result.size() == 3 ;
    if (passed) passed =  *result[0] == "a#" &&
		    *result[1] == "xxx!" && *result[2] == "space check";
    if (passed) 
	{
	    std::cout << "*** Illegal character function test: PASSED ***" 
		      << std::endl;
	}
    else
	{
	    std::cout << "*** Illegal character function test: FAILED ***" 
		      << std::endl;
	}
    

    std:: cout << std::endl;

    //---------------------------------------------------------------------------//

    // Test for acceptable lengths.

    std::cout << "String length utility test:" << std::endl;
    int low = 1;
    int high= 4;
    std::vector<VS_iter> result2 =
	rtt_dsxx::check_string_lengths(names.begin(), names.end(), low, high);
    if (result2.size() == 0) 
	std::cout << "All lengths OK!" << std::endl;
    else
    {
	std::cout << "*** Error in string definition -" << std::endl;
	for (size_t i=0; i<result2.size(); i++)
	    std::cout << "Size of string is not in allowable range: \"" 
		      << *result2[i] << "\"" << std::endl;
	std::cout << "Strings lengths must be greater than " << low 
		  << " and less than " << high << "." << std::endl;
    }

    passed = result2.size() == 2 ;
    if (passed) passed = *result2[0] == "space check" && *result2[1] == "";
    if (passed) 
	{
	    std::cout << "*** String length function test: PASSED ***" 
		      << std::endl;
	}
    else
	{
	    std::cout << "*** String length function test: FAILED ***" 
		      << std::endl;
	}

    std:: cout << std::endl;

    //---------------------------------------------------------------------------//

    // Test for unique names.

    std::cout << "Unique strings utility test:" << std::endl;
    std::vector<VS_iter> result3 =
	rtt_dsxx::check_strings_unique(names.begin(), names.end());
    if (result3.size() == 0) 
	std::cout << "All strings unique!" << std::endl;
    else
    {
	std::cout << "*** Error in string definition -" << std::endl;
	for (size_t i=0; i<result3.size(); i++)
	    std::cout << "Duplicate string found: \"" 
		      << *result3[i] << "\"" << std::endl;
	std::cout << "All strings must be unique!" << std::endl;
    }

    passed = result3.size() == 2 ;
    if (passed) passed = *result3[0] == "is" && *result3[1] == "test";
    if (passed) 
	{
	    std::cout << "*** Unique string function test: PASSED ***" 
		      << std::endl;
	}
    else
	{
	    std::cout << "*** Unique string function test: FAILED ***" 
		      << std::endl;
	}

    std:: cout << std::endl;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    //lint -e30 -e85 -e24 -e715 -e818 Suppress warnings about use of argv 
    //          (string comparison, unknown length, etc.)

    // version tag
    for( int arg = 1; arg < argc; arg++ )
	if( std::string(argv[arg]).find("--version") == 0 )
	{
	    std::cout << argv[0] << ": version " << rtt_dsxx::release() << std::endl; 
	    return 0;
	}

    try
    {
	// tests
	Check_Strings_test();
    }
    catch(rtt_dsxx::assertion const & ass)
    {
	std::cout << "Failure on Assertion: " << ass.what() << std::endl;
	return 1;
    }

    std::cout << "Done testing Check_Strings." << std::endl;

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstCheck_Strings.cc
//---------------------------------------------------------------------------//
