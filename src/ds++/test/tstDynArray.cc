//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstDynArray.cc
 * \author Kelly Thompson
 * \date   Fri Jan 30 17:53:51 2006
 * \brief  Unit tests and example usage for the DynArray class.
 * \note   Copyright 2006 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>

#include "../Assert.hh"
#include "../Release.hh"
#include "ds_test.hh"

#include "../DynArray.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_dynarray()
{
    // Create a DynArray using the default constructor.
    DynArray<int> myArray;

    // Add some data
    myArray[1]=100;
    myArray[2]=200;
    myArray[3]=300;
    myArray[15]=-1;

    if( myArray[1] != 100 ) ITFAILS;
    if( myArray[2] != 200 ) ITFAILS;
    if( myArray[3] != 300 ) ITFAILS;
    if( myArray[15] != -1 ) ITFAILS;

    cout << "myArray = " << myArray << endl;

    // Test Bounds
    {
        if( myArray.low() == 0 )
        {
            PASSMSG("Lower bound is zero.");
        }
        else
        {
            FAILMSG("Lower bound was not zero (the expected value).");
        }
        if( myArray.high() == 15 )
        {
            PASSMSG("Upper bound is 15.");
        }
        else
        {
            FAILMSG("Upper bound was not 15 (the expected value).");
        }
        int const expectedSize( 24 ); // provides for growth
        if( myArray.Get_size() == expectedSize )
        {
            ostringstream msg;
            msg << "The Get_size() member function returned the\n\t"
                << "expected value (" << expectedSize << ").";
            PASSMSG( msg.str() );
        }
        else
        {
            ostringstream msg;
            msg << "The Get_size() member function did not return the\n\t"
                << "expected value (" << expectedSize
                << "). Instead it provided the value:\n\t"
                << "myArray.Get_size() = " << myArray.Get_size();
            FAILMSG( msg.str() );
        }
    }
    
    // Copy Constructor
    {
        DynArray<int> copyOfArray( myArray );
        if( myArray == copyOfArray )
        {
            PASSMSG("Successful use of copy constructor.");
        }
        else
        {
            FAILMSG("Use of copy constructor failed.");
        }

        copyOfArray[-3] = 3;
        if( !(myArray == copyOfArray) )
        {
            PASSMSG("Successful use of operator==.");
        }
        else
        {
            FAILMSG("Use of operator== FAILS.");
        }        

        copyOfArray = myArray;
        copyOfArray[32];
        if( !(myArray == copyOfArray) )
        {
            PASSMSG("Successful use of operator==.");
        }
        else
        {
            FAILMSG("Use of operator== FAILS.");
        }        
    }

    // Assignment Operator
    {
        DynArray<int> secondArray;
        secondArray = myArray;
        if( myArray != secondArray )
        {
            FAILMSG("Use of assignment operator failed.");
        }
        else
        {
            PASSMSG("Successful use of assignment operator.");
        }

        // Manually set the bounds of what has been referenced.
        secondArray.low( 5 );
        secondArray.high( 5 );
        if( secondArray.low() != 5 ) ITFAILS;
        if( secondArray.high() != 5 ) ITFAILS;

        int dummyValue = secondArray[1];
        if( secondArray.low() != 1 ) ITFAILS;

        dummyValue = secondArray[13];
        if( secondArray.high() != 13 ) ITFAILS;
        
    }

    // Alternate constructor
    {
        int const size(5);
        int const base(1);
        int const init(99);
        float const growthFactor(1.5);
        DynArray<int> thirdArray( size, base, init, growthFactor );

        if( thirdArray[0] != init ) ITFAILS;
        if( thirdArray[3] != init ) ITFAILS;
        if( thirdArray.Get_defval() == init )
        {
            PASSMSG("Get_defval() member function returned the expected value.");
        }
        else
        {
            PASSMSG("Get_defval() member function did not return the expected value.");
        }
        if( thirdArray.Get_base() == base )
        {
            PASSMSG("Get_base() member function returned the expected value.");
        }
        else
        {
            PASSMSG("Get_base() member function did not return the expected value.");
        }
        if( thirdArray.Get_growthfactor() == growthFactor )
        {
            PASSMSG("Get_growthfactor() member function returned the expected value.");
        }
        else
        {
            PASSMSG("Get_growthfactor() member function did not return the expected value.");
        }

    }
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (std::string(argv[arg]) == "--version")
        {
            cout << argv[0] << ": version " 
                 << rtt_dsxx::release() << endl;
            return 0;
        }

    try
    {
        // >>> UNIT TESTS
        test_dynarray();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstDynArray, " 
                  << err.what()
                  << std::endl;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstDynArray, " 
		  << "An unknown exception was thrown"
                  << std::endl;
        return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
              << std::endl;
    if (rtt_ds_test::passed) 
    {
        std::cout << "**** tstDynArray Test: PASSED"
                  << std::endl;
    }
    std::cout <<     "*********************************************" 
              << std::endl;
    std::cout << std::endl;
    

    std::cout << "Done testing tstDynArray"
              << std::endl;
    
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstDynArray.cc
//---------------------------------------------------------------------------//
