//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/tDummyOpacity.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  9 15:50:53 2001
 * \brief  GrayOpacity and Multigroup opacity test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_test.hh"
#include "DummyGrayOpacity.hh"
#include "../Release.hh"
#include "../GrayOpacity.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

using rtt_cdi_test::match;
using rtt_cdi::GrayOpacity;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void gray_opacity_test()
{
    // ---------------------------- //
    // Create a GrayOpacity object. //
    // ---------------------------- //

    SP< GrayOpacity > spDGO;
    
    if ( spDGO = new rtt_cdi_test::DummyGrayOpacity() )
    {
	PASSMSG("SP to new GrayOpacity object created.");
    }
    else
    {
	FAILMSG("Unable to create a SP to new GrayOpacity object.");
    }

    // ------------------------ //
    // Dummy Gray Opacity Tests //
    // ------------------------ //
    
    double temperature          = 0.1;                          // keV
    double density              = 27.0;                         // g/cm^3
    double tabulatedGrayOpacity = temperature + density/1000.0; // cm^2/g
    
    double opacity = spDGO->getOpacity( temperature, density );

    if ( match ( opacity, tabulatedGrayOpacity ) ) 
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity computation was good.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity value is out of spec.";
	FAILMSG( message.str() );
    }

    // try using a vector of temps.

    std::vector< double > vtemperature(2);
    vtemperature[0] = 0.5;  // keV
    vtemperature[1] = 0.7;  // keV
    density         = 0.35; // g/cm^3
    std::vector< double > vRefOpacity( vtemperature.size() );
    for ( int i=0; i<vtemperature.size(); ++i )
	vRefOpacity[i] = vtemperature[i] + density/1000;

    std::vector< double > vOpacity = spDGO->getOpacity( vtemperature, 
							density );  

    if ( match ( vOpacity, vRefOpacity ) ) 
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity computation was good for a vector of temps.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity value is out of spec. for a vector of temps.";
	FAILMSG( message.str() );
    }
	

    // try using a vector of densities.

    std::vector< double > vdensity(5);
    vdensity[0] = 0.5;
    vdensity[1] = 1.0;
    vdensity[2] = 3.3;
    vdensity[3] = 5.0;
    vdensity[4] = 27.0;

    vRefOpacity.resize( vdensity.size() );
    for ( int i=0; i<vdensity.size(); ++i )
	vRefOpacity[i] = temperature + vdensity[i]/1000;

    vOpacity = spDGO->getOpacity( temperature, vdensity );

    if ( match ( vOpacity, vRefOpacity ) )
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity computation was good"
		<< " for a vector of densities.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity value is out of spec."
		<< " for a vector of densities.";
	FAILMSG( message.str() );
    }
}

//---------------------------------------------------------------------------//

void multigroup_opacity_test()
{
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_cdi::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	gray_opacity_test();
	multigroup_opacity_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tDummyOpacity, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_test::passed) 
    {
        cout << "**** tDummyOpacity Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tDummyOpacity." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tDummyOpacity.cc
//---------------------------------------------------------------------------//
