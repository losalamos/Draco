//----------------------------------*-C++-*----------------------------------//
// tQuadrature.cc
// Kelly Thompson
// Thu Mar 2 13:07:00 2000
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// revision history:
// -----------------
// 1.1) Original
// ...  ...
// 1.5) Added a pass() call for each quadrature set tested.
//      QuadCreator.quadCreate member function name changed to lower
//         case.
//      Moved "using" statements inside of namespace.
//      Added "using std::fabs".

#include "tQuadrature.hh"
#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"

#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	using rtt_quadrature_test::tQuadrature;
	return rtt_dsxx::SP<TestApp> ( new tQuadrature( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame


// tQuadrature Stuff
//----------------------------------------
namespace rtt_quadrature_test {

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::fabs;
using rtt_dsxx::SP;
using rtt_quadrature::QuadCreator;
using rtt_quadrature::Quadrature;

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
    
    // we will only look at S4 Sets in this test.
    const int sn_order = 4;

    // total number of quadrature sets to be tested.
    const int nquads = 3;

    // Declare an enumeration object that specifies the Quadrature set to be
    // tested.

    // Quadrature sets to be tested:
    //
    // #   Qid        Description
    // -   --------   ------------
    // 0   GaussLeg   1D Gauss Legendre
    // 1   LevelSym2D 2D Level Symmetric
    // 2   LevelSym   3D Level Symmetric

    QuadCreator::Qid qid[nquads] = { QuadCreator::GaussLeg,
				     QuadCreator::LevelSym2D,
				     QuadCreator::LevelSym };

    // mu0 holds mu for the first direction for each quadrature set tested.
    double mu0[nquads] = { 0.8611363116,
			  -0.350021174581541, -0.350021174581541 };
    
    SP<Quadrature> spQuad;

    // loop over quadrature types to be tested.

    for ( int ix = 0; ix < nquads; ++ix ) {
	
	// Verify that the enumeration value matches its int value.
	if ( qid[ix] != ix ) {
	    fail() << "Setting QuadCreator::Qid enumeration failed.";
	    break;
	} else {
	    // Instantiate the quadrature object.
	    spQuad = QuadratureCreator.quadCreate( qid[ix], sn_order ); 

	    // print the name of the quadrature set that we are testing.
	    string qname = spQuad->name();
	    cout << "\nTesting the "  << qname
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
		// get the total omega object.
		vector< vector<double> > omega = spQuad->getOmega();
		// compare values.
		if ( mu.size() != spQuad->getNumAngles() )
		    fail() << "The direction vector has the wrong length.";
		else if ( fabs( mu[0] + mu0[ix] ) >= TOL ) 
		    fail() << "mu[0] has the wrong value."; 
		else if ( fabs( mu[1] - omega_1[0] ) >= TOL )
		    fail() << "mu[1] != omega_1[0].";
		else if ( fabs( omega[1][0] - omega_1[0] ) >= TOL )
		    fail() << "omega[1][0] != omega_1[0].";
		else {
		    spQuad->display();
		    cout << endl << endl; // end of this quadrature type
		}
	    }
	    pass() << "Passed all tests for the " << qname << " quadrature set.";
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

