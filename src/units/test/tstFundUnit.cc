//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   units/test/tstFundUnit.cc
 * \author Kelly Thompson
 * \date   Wed Oct  8 13:50:19 2003
 * \brief  
 * \note   Copyright (C) 2003-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>
#include <sstream>

#include "ds++/Soft_Equivalence.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "../FundUnit.hh"
#include "units_test.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

void test_construction()
{
    using std::cout;
    using std::endl;
    using std::string;
    using std::ostringstream;
    using rtt_dsxx::soft_equiv;
    using rtt_units::FundUnit;
    using rtt_units::Ltype;
    using rtt_units::L_null;
    using rtt_units::L_labels;
    using rtt_units::L_cf;

    FundUnit< Ltype > myLength( L_null, L_cf, L_labels );

    // Is the enum value correct?
    if( myLength.enumVal() == L_null )
    {
	PASSMSG("Enum Value was set correctly.");
    }
    else
    {
	ostringstream msg;
	msg << "Enum Value was not set correctly "
	    << "(" << myLength.enumVal() << " != L_null (0) )" << endl;
	FAILMSG(msg.str());
    }

    // Is the conversion factor correct?
    if( soft_equiv( myLength.cf(), 0.0 ) )
    {
	PASSMSG("Conversion Factor was set correctly.");
    }
    else
    {
	ostringstream msg;
	msg << "Conversion Factor was not set correctly "
	    << "(" << myLength.cf() << " != 0.0)" << endl;
	FAILMSG(msg.str());
    }

    // Is the label correct?
    if( myLength.label() == string( "NA" ) )
    {
	PASSMSG("Unit label was set correctly.");
    }
    else
    {
	ostringstream msg;
	msg << "Unit label was not set correctly "
	    << "(\"" << myLength.label() << "\" != \"NA\")" << endl;
	FAILMSG(msg.str());
    }

    //-----------------------------------//
    
    // Use "long labels" string...

    FundUnit< Ltype > myLength2( rtt_units::L_cm, 
				 rtt_units::L_cf, 
				 rtt_units::L_long_labels );

    // Is the label correct?
    std::string sentinelValue( "centimeter" );
    if( myLength2.label() == sentinelValue )
    {
	PASSMSG("Unit long label was set correctly.");
    }
    else
    {
	ostringstream msg;
	msg << "Unit long label was not set correctly "
	    << "(\"" << myLength2.label() << "\" != \"" << sentinelValue 
	    << "\")" << endl;
	FAILMSG(msg.str());
    }
    // Is the conversion factor correct?
    if( soft_equiv( myLength2.cf(), 100.0 ) )
    {
	PASSMSG("Conversion Factor was set correctly.");
    }
    else
    {
	ostringstream msg;
	msg << "Conversion Factor was not set correctly "
	    << "(" << myLength2.cf() << " != 100.0)" << endl;
	FAILMSG(msg.str());
    }

    return;

}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;
    using std::string;
    using rtt_dsxx::release;

    // version tag
    for( int arg = 1; arg < argc; arg++ )
    {
	if( string( argv[arg] ) == "--version" )
	{
	    cout << argv[0] << ": version " << release() << endl;
	    return 0;
	}
    }

    try  //  >>>>> UNIT TESTS <<<<<
    {
	test_construction();
    }
    catch( rtt_dsxx::assertion &assert )
    {
	cout << "While testing tstFundUnit, " << assert.what() << endl;
	return 1;
    }

    // status of test
    cout << endl << "*********************************************" << endl;
    if( rtt_units_test::passed )
        cout << "**** tstFundUnit Test: PASSED" << endl;
    else
        cout << "**** tstFundUnit Test: FAILED" << endl;
    cout <<         "*********************************************" << endl << endl;
    
    cout << "Done testing tstFundUnit." << endl;
    return 0;
}   

//---------------------------------------------------------------------------//
// end of tstFundUnit.cc
//---------------------------------------------------------------------------//
