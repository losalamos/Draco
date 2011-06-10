//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstData_Table.cc
 * \author Paul Henning
 * \brief  DBC_Ptr tests.
 * \note   Copyright (c) 1997-2011 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds_test.hh"
#include "../Release.hh"
#include "../Data_Table.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <typeinfo>
#include <sstream>

using namespace std;
using rtt_dsxx::Data_Table;


//---------------------------------------------------------------------------//

void test_array()
{
    const int array[3] = {10,11,12};

    bool caught = false;
    try
    {
        Data_Table<int> dt(array, array+3);
    }
    catch(rtt_dsxx::assertion &/* error */)
    {
        caught = true;
    }
    if(caught) ITFAILS;

    caught = false;
    try
    {
        Data_Table<int> dt(array+3,array);
    }
    catch(rtt_dsxx::assertion & /* error */)
    {
        caught = true;
    }
    if(!caught) ITFAILS;


    Data_Table<int> dt(array, array+3);

    caught = false;
    try
    {
        if(dt.size() != 3) ITFAILS;
        if(dt[0] != array[0]) ITFAILS;
        if(dt[1] != array[1]) ITFAILS;
        if(dt[2] != array[2]) ITFAILS;
        if(dt.front() != array[0]) ITFAILS;
        if(dt.back() != array[2]) ITFAILS;
        
        if(dt.begin() != array) ITFAILS;
        if(dt.end() != array + 3) ITFAILS;

        if (dt.access() != &dt[0]) ITFAILS;

        {
            Data_Table<int> dt3(dt);
            if(dt3.size() != dt.size()) ITFAILS;
            if(dt3.begin() != dt.begin()) ITFAILS;
            if(dt3.end() != dt.end()) ITFAILS;
            if(dt3.front() != dt.front()) ITFAILS;
            if(dt3.back() != dt.back()) ITFAILS;
        }


        {
            Data_Table<int> dt3;
            if(dt3.size() != 0) ITFAILS;
            dt3 = dt;
            if(dt3.size() != dt.size()) ITFAILS;
            if(dt3.begin() != dt.begin()) ITFAILS;
            if(dt3.end() != dt.end()) ITFAILS;
            if(dt3.front() != dt.front()) ITFAILS;
            if(dt3.back() != dt.back()) ITFAILS;
        }
    }
    catch(rtt_dsxx::assertion & /* error */) 
    {
        caught = true;
    }
    if(caught) ITFAILS;


#ifdef DEBUG
/*
  GCC will issue a warning at compile time for a Release build (with
  -ftree-vrp, which is enabled by default with -O2 or higher).  The warning
  appears because the size of dt is known at compile time.
*/
    caught = false;
    try
    {
        std::cout << dt[3];
    }
    catch(rtt_dsxx::assertion & /* error */) 
    {
        caught = true;
    }
    if(!caught) ITFAILS;
#endif

    if (rtt_ds_test::passed)
        PASSMSG("test_array");
    else
        FAILMSG("test_array FAILED!");
}

//---------------------------------------------------------------------------//

void test_scalar()
{
    Data_Table<int> dt(32);

    bool caught = false;
    try
    {
        if(dt.size() != 1) ITFAILS;
        if(dt[0] != 32) ITFAILS;

        Data_Table<int> dt2;
        dt2 = dt;
        if(dt[0] != dt2[0]) ITFAILS;
        if(&(dt[0]) == &(dt2[0])) ITFAILS;
        if(dt.front() != 32) ITFAILS;
        if(dt.back() != 32) ITFAILS;
        if(*(dt.begin()) != 32) ITFAILS;
        if((dt.end() - dt.begin()) != 1) ITFAILS;

        Data_Table<int> dt3(dt2);
        if(dt[0] != dt3[0]) ITFAILS;
        if(&(dt[0]) == &(dt3[0])) ITFAILS;

        dt3 = dt3;
        if(dt[0] != dt3[0]) ITFAILS;
        if(&(dt[0]) == &(dt3[0])) ITFAILS;
    }
    catch(rtt_dsxx::assertion & /* error */) 
    {
        caught = true;
    }
    if(caught) ITFAILS;


    caught = false;
    try
    {
        std::cout << dt[1];
    }
    catch(rtt_dsxx::assertion & /* error */) 
    {
        caught = true;
    }
    if(!caught) ITFAILS;

    if (rtt_ds_test::passed)
        PASSMSG("test_scalar");
    else
        FAILMSG("test_scalar FAILED!");
}


//---------------------------------------------------------------------------//
/*!
  As noted in the Class Declaration/Definition, this use of Data_Table is no
  longer supported because it did not conform to the C++ standard.
*/
  
// void test_vector()
// {
//     std::vector<char> v;
//     v.push_back('a');
//     v.push_back('b');
//     v.push_back('c');

//     Data_Table<char> dt(v);

//     bool caught = false;
//     try
//     {
//         if(dt.size() != 3) ITFAILS;
//         if(dt[0] != 'a') ITFAILS;
//         if(dt[1] != 'b') ITFAILS;
//         if(dt[2] != 'c') ITFAILS;

//         Data_Table<char> dt2;
//         dt2 = dt;
//         if(dt[0] != dt2[0]) ITFAILS;
//         if(&(dt[0]) != &(dt2[0])) ITFAILS;
//         if(dt.front() != 'a') ITFAILS;
//         if(dt.back() != 'c') ITFAILS;
//         if(*(dt.begin()) != 'a') ITFAILS;
//         if((dt.end() - dt.begin()) != 3) ITFAILS;
//     }
//     catch(rtt_dsxx::assertion &ass) 
//     {
//         caught = true;
//     }
//     if(caught) ITFAILS;


//     caught = false;
//     try
//     {
//         std::cout << dt[3];
//     }
//     catch(rtt_dsxx::assertion &ass) 
//     {
//         caught = true;
//     }
//     if(!caught) ITFAILS;

//     if (rtt_ds_test::passed)
// 	PASSMSG("test_vector");
//     else
// 	FAILMSG("test_vector FAILED!");
// }

//---------------------------------------------------------------------------//

int 
main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() 
		 << endl;
	    return 0;
	}

#if DBC
    try
    {
	// >>> UNIT TESTS
	test_array();
        test_scalar();
//        test_vector();
    }
    catch (rtt_dsxx::assertion &error)
    {
        cout << "While testing tstData_Table_Ptr, " << error.what()
             << endl;
        return 1;
    }

    catch (...)
    {
        cout << "caught uncaught exception" << std::endl;
        return 10;
    }

#endif

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_ds_test::passed) 
    {
        cout << "**** tstData_Table Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstData_Table." << endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstData_Table.cc
//---------------------------------------------------------------------------//
