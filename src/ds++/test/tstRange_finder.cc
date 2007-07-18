//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   uncleMcflux/test/tstRange_finder.cc
 * \author Mike Buksas
 * \date   Thu Feb  6 12:43:22 2003
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds_test.hh"
#include "../Release.hh"
#include "../Assert.hh"
#include "../Range_Finder.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <vector>

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_range_finder_left()
{

    double v[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    vector<double> values(v,v+10);

    int index = Range_finder_left(v,v+10, 1.5);   
    if (index != 1) ITFAILS;

    index = Range_finder_left(values.begin(), values.end(), 2.5);
    if (index != 2) ITFAILS;

    // Check for equality at all values:
    for (int i=0; i < 10; ++i)
    {
	index = Range_finder_left(v, v+10, static_cast<double>(i));
	if (index != i) ITFAILS;
    }

    // For equality with the last value, we should get n-1 with end catching:
    index = Range_finder_left_catch_end(v,v+10, 9.0);
    if (index != 8) ITFAILS;

//     index = Range_finder_left(v,v+10, 42.69);
//     if (index != -1) ITFAILS;

//     index = Range_finder_left(v+5,v+10, 1.0);
//     if (index != -1) ITFAILS;

    double rv[10] = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0};

    vector<double> rvalues(rv,rv+10);

    index = Range_finder(rvalues.rbegin(), rvalues.rend(), 5.5, LEFT);
    if (index != 5) ITFAILS;

    index = Range_finder_left(rvalues.rbegin(), rvalues.rend(), 5.0);
    if (index != 5) ITFAILS;

//     index = Range_finder_left(rvalues.rbegin(), rvalues.rend(), 10.12);
//     if (index != -1) ITFAILS;

   
    if (rtt_ds_test::passed)
	PASSMSG("tstRange_finder_left PASSED");

}

void test_range_finder_right()
{

    double v[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    int index;

    // Check for equality at all values. Note that 0 comes back as interval
    // -1 (e.g. out of range).
    for (int i=1; i < 10; ++i)
    {
	index = Range_finder(v, v+10, static_cast<double>(i), RIGHT);
	if (index != i-1) ITFAILS;
    }

    index = Range_finder_right_catch_end(v, v+10, 0.0);
    if (index != 0) ITFAILS;

    if (rtt_ds_test::passed)
	PASSMSG("tstRange_finder_right PASSED");
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	test_range_finder_left();

	test_range_finder_right();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstRange_finder, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_ds_test::passed) 
    {
        cout << "**** tstRange_finder Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstRange_finder." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstRange_finder.cc
//---------------------------------------------------------------------------//
