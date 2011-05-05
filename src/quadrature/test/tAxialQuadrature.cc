//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tAxialQuadrature.cc
 * \author Jae Chang
 * \date   Fri Oct 08 10:26:41 2004
 * \brief  quadrature package test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "ds++/Release.hh"

#include "quadrature_test.hh"

using namespace std;

using rtt_quadrature::QuadCreator;
using rtt_quadrature::Quadrature;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
/*!
 * \brief Tests the Quadcrator and Quadtrature constructors and access
 * routines. 
 *
 * To add a quadrature to this test the following items must be changed: add
 * new enumeration to Qid[] array.  Verify nquads is set to the correct number
 * of quadrature sets being tested.
 */
void quadrature_test()
{
    using std::ostringstream;
    using std::endl;
    using std::cout;

    // create an object that is responsible for creating quadrature objects.
    QuadCreator QuadratureCreator;
    
    // we will only look at S2 in this test.
    const int sn_order = 2;

    // total number of quadrature sets to be tested.
    const int nquads = 1;

    // Declare an enumeration object that specifies the Quadrature set to be
    // tested.

    // Quadrature sets to be tested:
    //
    // #   Qid        Description
    // -   --------   ------------
    // 0   Axial1D    1D Axial

    QuadCreator::Qid qid[nquads] = { QuadCreator::Axial1D };

    SP< const Quadrature > spQuad;

    // loop over quadrature types to be tested.

    for ( int ix = 0; ix < nquads; ++ix ) {
	
	// Verify that the enumeration value matches its int value.
	if ( qid[ix] != QuadCreator::Axial1D ) {
	    FAILMSG("Setting QuadCreator::Qid enumeration failed.");
	    break;
	} else {
	    // Instantiate the quadrature object.
	    spQuad = QuadratureCreator.quadCreate( qid[ix], sn_order ); 

	    // print the name of the quadrature set that we are testing.
	    string qname = spQuad->name();
	    cout << "\nTesting the "  << qname
		 << "Quadrature set." << endl;
	    cout << "   Sn Order         = " << spQuad->getSnOrder() << endl;
	    cout << "   Number of Ordinates = " << spQuad->getNumOrdinates() << endl;
            cout << "   Parser Name = " << spQuad->parse_name() << endl;
            cout << "   Class = " << spQuad->getClass() << endl;

	    // If the object was constructed sucessfully then we continue
	    // with the tests.
	    if ( ! spQuad )
		FAILMSG("QuadCreator failed to create a new quadrature set.")
	    else {
		// get the mu vector
		vector<double> mu = spQuad->getMu();
		if ( mu.size() != spQuad->getNumOrdinates() )
		    FAILMSG("The direction vector has the wrong length.")
		else 
		{
		    spQuad->display();
		    cout << endl << endl; // end of this quadrature type
		}
	    }
	    std::ostringstream msg;
	    msg << "Passed all tests for the " << qname 
		<< " quadrature set.";
	    PASSMSG( msg.str() );
	}
    }
    
    // Test dimensionality accessor.
    size_t const expected_dim( 1 );
    if( spQuad->dimensionality() == expected_dim )
    {
	PASSMSG("Found expected dimensionality value.");
    }
    else
    {
	ostringstream msg;
	msg << "Did not find expected dimensionality == " << expected_dim << "."
	    << endl
	    << "spQuad returned dim=" << spQuad->dimensionality() << " instead."
	    << endl;	
	FAILMSG(msg.str()); 
    }

    return;
} // end of quadrature_test

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
	quadrature_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tAxialQuadrature, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_quadrature_test::passed) 
    {
        cout << "**** tAxialQuadrature Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tAxialQuadrature." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tQuadrature.cc
//---------------------------------------------------------------------------//
