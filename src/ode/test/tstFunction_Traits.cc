//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ode/test/tstFunction_Traits.cc
 * \author Kent Budge
 * \date   Wed Aug 18 10:28:16 2004
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <typeinfo>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "../Function_Traits.hh"

using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

class Test_Functor
{
  public:
    typedef double return_type;
};

void tstFunction_Traits( UnitTest & ut )
{
    if (typeid(Function_Traits<double (*)(double)>::return_type) != 
	typeid(double))
    {
	ut.failure("return_type NOT correct");
    }
    else
    {
	ut.passes("return_type correct");
    }

    if (typeid(Function_Traits<Test_Functor>::return_type) != 
	typeid(Test_Functor::return_type))
    {
	ut.failure("return_type NOT correct");
    }
    else
    {
	ut.passes("return_type correct");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        ScalarUnitTest ut( argc, argv, release );
	tstFunction_Traits(ut);
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown. " << endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstFunction_Traits.cc
//---------------------------------------------------------------------------//
