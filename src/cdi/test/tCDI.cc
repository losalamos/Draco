//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/tCDI.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  9 15:52:01 2001
 * \brief  CDI test executable.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_test.hh"
#include "DummyGrayOpacity.hh"
#include "DummyMultigroupOpacity.hh"
#include "DummyEoS.hh"
#include "../Release.hh"
#include "../CDI.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <typeinfo>

using namespace std;

using rtt_cdi_test::match;
using rtt_cdi_test::DummyGrayOpacity;
using rtt_cdi_test::DummyMultigroupOpacity;
using rtt_cdi_test::DummyEoS;
using rtt_cdi::CDI;
using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_cdi::EoS;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void check_CDI(const CDI &cdi)
{
    // check cdi, note that the different combinations of rtt_cdi::Model and
    // rtt_cdi::Reaction will yield the same results because DummyOpacity,
    // DummyMultigroupOpacity, and DummyEoS all yield the same stuff.  These
    // have all been tested in tDummyOpacity and tDummyEoS.  Here we just
    // check the types
    
    if (typeid(*cdi.gray(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)) == 
	typeid(DummyGrayOpacity))
    {
	PASSMSG("CDI gray() returned the correct type!");
    }
    else
    {
	FAILMSG("CDI gray() did not return the correct type!");
    }
    
    if (typeid(*cdi.gray(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)) == 
	typeid(DummyGrayOpacity))
    {
	PASSMSG("CDI gray() returned the correct type!");
    }
    else
    {
	FAILMSG("CDI gray() did not return the correct type!");
    }
    
    if (typeid(*cdi.mg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)) == 
	typeid(DummyMultigroupOpacity))
    {
	PASSMSG("CDI mg() returned the correct type!");
    }
    else
    {
	FAILMSG("CDI mg() did not return the correct type!");
    }

    
    if (typeid(*cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)) == 
	typeid(DummyMultigroupOpacity))
    {
	PASSMSG("CDI mg() returned the correct type!");
    }
    else
    {
	FAILMSG("CDI mg() did not return the correct type!");
    }
    
    if (typeid(*cdi.eos()) == typeid(DummyEoS))
    {
	PASSMSG("CDI eos() returned the correct type!");
    }
    else
    {
	FAILMSG("CDI eos() did not return the correct type!");
    }

    // gray test case: Find the value of opacity at T=0.35 keV and rho = 27.2
    // g/cm^3.  For DummyGrayOpacity the value should be .35272 cm^2/g.
     
    double temp       = 0.35;               // keV
    double dens       = 27.2;               // g/cm^3
    double refOpacity = temp + dens/1000.0; // cm^2/g

    double opacity = cdi.gray(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION)->
	getOpacity(temp, dens);

    if (match(opacity, refOpacity))
    {
	PASSMSG("CDI.gray()->getOpacity is ok.");
    }
    else
    {
	FAILMSG("CDI.gray()->getOpacity is not ok.");
    }

    // mg test case: Find the mg opacities at T=0.35 keV and rho = 27.2
    // g/cm^3.  For DummyMultigroupOpacity the values should be { }.  Three
    // groups    
    int ng = 3;

    // The energy groups in DummyMultigroupOpacity are hardwired to
    // be { 0.05, 0.5, 5.0, 50.0 } keV.
    std::vector< double > energyBoundary(4);
    energyBoundary[0] = 0.05;
    energyBoundary[1] = 0.5;
    energyBoundary[2] = 5.0;
    energyBoundary[3] = 50.0;

    std::vector< double > vRefOpacity( 3 );
    for ( int ig=0; ig<ng; ++ig )
	vRefOpacity[ig] = 2.0*(temp+dens/1000.0)
	    /(energyBoundary[ig]+energyBoundary[ig+1]);

    std::vector< double > vOpacity( ng );
    vOpacity = cdi.mg(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING)->
	getOpacity( temp, dens );

    if ( match( vOpacity, vRefOpacity ) )
    {
	PASSMSG("CDI->mg()->getOpacity(T,rho) is ok.");
    }
    else
    {
	FAILMSG("CDI->mg()->getOpacity(T,rho) is not ok.");
    }

    // Test the EoS plug-in component of this CDI.

    double refCve = temp + dens/1000.0;
	    
    double Cve = cdi.eos()->getElectronHeatCapacity( temp, dens );

    if ( match( Cve, refCve ) )
    {
	ostringstream message;
	message << "CDI->eos()->getElectronHeatCapacity( dbl, dbl )\n\t"
		<< "returned a value that matched the reference value.";
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "CDI->eos()->getElectronHeatCapacity( dbl, dbl )\n\t"
		<< "returned a value that was out of spec.";
	FAILMSG(message.str());
    }
}

//---------------------------------------------------------------------------//

void test_CDI()
{
    // make SPs to opacity and EoS objects
    SP<const GrayOpacity>       gray_planck_abs;
    SP<const GrayOpacity>       gray_iso_scatter;
    SP<const MultigroupOpacity> mg_planck_abs;
    SP<const MultigroupOpacity> mg_iso_scatter;
    SP<const EoS>               eos;

    // assign to dummy state objects
    gray_planck_abs  = new DummyGrayOpacity(rtt_cdi::ABSORPTION,
					    rtt_cdi::PLANCK);
    gray_iso_scatter = new DummyGrayOpacity(rtt_cdi::SCATTERING,
					    rtt_cdi::ISOTROPIC);

    mg_planck_abs    = new DummyMultigroupOpacity(rtt_cdi::ABSORPTION,
						  rtt_cdi::PLANCK);
    mg_iso_scatter   = new DummyMultigroupOpacity(rtt_cdi::SCATTERING,
						  rtt_cdi::ISOTROPIC);
    
    eos              = new DummyEoS();
    
    // make a CDI, it should be empty
    CDI cdi;
    for (int i = 0; i < rtt_cdi::constants::num_Models; i++)
	for (int j = 0; j < rtt_cdi::constants::num_Reactions; j++)
	{
	    // casting these is inherently dangerous, but because this is a
	    // test I won't go nuts
	    rtt_cdi::Model m    = static_cast<rtt_cdi::Model>(i);
	    rtt_cdi::Reaction r = static_cast<rtt_cdi::Reaction>(j);

	    if (cdi.isGrayOpacitySet(m,r))       ITFAILS;
	    if (cdi.isMultigroupOpacitySet(m,r)) ITFAILS;
	}
    if (cdi.isEoSSet()) ITFAILS;

    // now assign stuff to it
    cdi.setGrayOpacity(gray_planck_abs);
    cdi.setGrayOpacity(gray_iso_scatter);
    cdi.setMultigroupOpacity(mg_planck_abs);
    cdi.setMultigroupOpacity(mg_iso_scatter);
    cdi.setEoS(eos);

    // make sure these are assigned
    if (cdi.isGrayOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
	PASSMSG("Gray planck absorption set!");
    }
    else
    {
	FAILMSG("Gray planck absorption not set!");
    }

    if (cdi.isGrayOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
	PASSMSG("Gray isotropic scattering set!");
    }
    else
    {
	FAILMSG("Gray isotropic scattering not set!");
    }

    if (cdi.isMultigroupOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
	PASSMSG("Multigroup planck (in-group) absorption set!");
    }
    else
    {
	FAILMSG("Multigroup planck (in-group) absorption not set!");
    }

    if (cdi.isMultigroupOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
	PASSMSG("Multigroup isotropic scattering set!");
    }
    else
    {
	FAILMSG("Multigroup isotropic scattering not set!");
    }

    if (cdi.isEoSSet())
    {
	PASSMSG("EoS set!");
    }
    else
    {
	FAILMSG("EoS not set!");
    }

    // catch some exceptions
    bool caught = false;
    try 
    {
	cdi.setGrayOpacity(gray_planck_abs);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	PASSMSG("Good, we caught an overwrite exception.");
	caught = true;
    }
    if (!caught) 
    {
	FAILMSG("Failed to catch overset exception!");
    }

    caught = false;
    try
    {
	cdi.mg(rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	PASSMSG("Good, we caught an illegal access exception.");
	caught = true;
    }
    if (!caught)
    {
	FAILMSG("Failed to catch an illegal access exception!");
    }

    // check the cdi through a function call
    check_CDI(cdi);

    // reset and make sure we are empty
    cdi.reset();
    if (!cdi.isGrayOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
	PASSMSG("Gray planck absorption unset!");
    }
    else
    {
	FAILMSG("Gray planck absorption is still set!");
    }

    if (!cdi.isGrayOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
	PASSMSG("Gray isotropic scattering unset!");
    }
    else
    {
	FAILMSG("Gray isotropic scattering is still set!");
    }

    if (!cdi.isMultigroupOpacitySet(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION))
    {
	PASSMSG("Multigroup planck (in-group) absorption unset!");
    }
    else
    {
	FAILMSG("Multigroup planck (in-group) absorption is still set!");
    }

    if (!cdi.isMultigroupOpacitySet(rtt_cdi::ISOTROPIC, rtt_cdi::SCATTERING))
    {
	PASSMSG("Multigroup isotropic scattering unset!");
    }
    else
    {
	FAILMSG("Multigroup isotropic scattering is still set!");
    }

    if (!cdi.isEoSSet())
    {
	PASSMSG("EoS unset!");
    }
    else
    {
	FAILMSG("EoS is still set!");
    }

    caught = false;
    try
    {
	cdi.mg(rtt_cdi::PLANCK, rtt_cdi::ABSORPTION);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	PASSMSG("Good, we caught an illegal access exception.");
	caught = true;
    }
    if (!caught)
    {
	FAILMSG("Failed to catch an illegal access exception!");
    }
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
	test_CDI();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tCDI, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_test::passed) 
    {
        cout << "**** tCDI Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tCDI." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tCDI.cc
//---------------------------------------------------------------------------//
