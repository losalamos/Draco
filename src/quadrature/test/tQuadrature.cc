//----------------------------------*-C++-*----------------------------------//
// tQuadrature.cc
// Kelly Thompson
// Thu Mar 2 13:07:00 2000
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "tQuadrature.hh"
#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"

#include <iostream>
#include <vector>
#include <sstream>

using rtt_dsxx::SP;
using rtt_quadrature::QuadCreator;
using rtt_quadrature::Quadrature;

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame {
    SP<TestApp> TestApp::create( int &argc, char *argv[],
				 std::ostream& os_in ) {
	    using rtt_quadrature_test::tQuadrature;
	    return SP<TestApp> ( new tQuadrature( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame


// tQuadrature Stuff
//----------------------------------------

namespace rtt_quadrature_test {

using std::vector;
using std::string;
using std::cout;
using std::endl;

tQuadrature::tQuadrature( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
{
    os() << "Created tQuadrature" << endl;
}

string tQuadrature::version() const
{
    return rtt_quadrature::release();
}


// To add a quadrature to this test the following items must be changed:
//    add new enumeration to Qid[] array.
//    add new mu[0] value to mu0[] array.
//    verify nquads is set to the correct number of quadrature sets being
//       tested. 
string tQuadrature::runTest()
{
    // double precesion values will be tested for correctness against this
    // tolerance. 
    const double TOL = 1.0e-10; 

    // create an object that is responsible for creating quadrature objects.
    QuadCreator QuadratureCreator;
//     pass(" Instantiate QuadCreator") 
// 	<< "Successfully instantiated a QuadCreator object." << endl;
    
    // we will only look at S4 Sets in this test.
    int sn_order = 4;

    // total number of quadrature sets to be tested.
    const int nquads = 2;

    // Declare an enumeration object that specifies the Quadrature set to be
    // tested.

    // Quadrature sets to be tested:
    //
    // #   Qid        Description
    // -   --------   ------------
    // 0   GaussLeg   1D Gauss Legendre
    // 1   LevelSym   3D Level Symmetric

    QuadCreator::Qid qid[nquads] = { QuadCreator::GaussLeg,
				     QuadCreator::LevelSym  };

    // mu0 holds mu for the first direction for each quadrature set tested.
    double mu0[nquads] = { 0.8611363116,
			  -0.350021174581541 };
    
    SP<Quadrature> spQuad;

    // loop over quadrature types to be tested.

    for ( int ix = 0; ix < nquads; ++ix ) {
	
	// Verify that the enumeration value matches its int value.
	if ( qid[ix] != ix ) {
	    fail() << "Setting QuadCreator::Qid enumeration failed.";
	    break;
	} else {
	    // Instantiate the quadrature object.
	    spQuad = QuadratureCreator.QuadCreate( qid[ix], sn_order ); 

	    // print the name of the quadrature set that we are testing.
	    cout << "\nTesting the "  << spQuad->name() 
		 << "Quadrature set." << endl;
	    cout << "   Sn Order         = " << spQuad->getSnOrder() << endl;
	    cout << "   Number of Angles = " << spQuad->getNumAngles() << endl;

	    // If the object was constructed sucessfully then we continue
	    // with the tests.
	    if ( ! spQuad )
		fail() << "QuadCreator failed to create a new quadrature set.";
	    else {
		// get the mu vector
		vector<double> mu = spQuad->getMu();
		// get the omega vector for direction m=1.
		vector<double> omega_1 = spQuad->getOmega(1);
		if ( mu.size() != spQuad->getNumAngles() )
		    fail() << "The direction vector has the wrong length.";
		else if ( fabs( mu[0] + mu0[ix] ) >= TOL ) 
		    fail() << "mu[0] has the wrong value."; 
		else if ( fabs( mu[1] - omega_1[0] ) >= TOL )
		    fail() << "mu[1] != omega_1[0].";
		else {
		    spQuad->display();
		    cout << endl << endl; // end of this quadrature type
		}
	    }
	}
    }

    // Print the test result.
    // ----------------------------------------

    if (passed()) {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}

} // end namespace rtt_Quadrature_test


//---------------------------------------------------------------------------//
//                            end of tQuadrature.cc
//---------------------------------------------------------------------------//

