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
using std::cout;
using std::endl;
#include <vector>
#include <sstream>

static double TOL = 1.0e-14;


// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame {
// what does this do?
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
				       std::ostream& os_in ) 
	{
	    using rtt_dsxx::SP;
	    using rtt_quadrature_test::tQuadrature;
	    
	    return SP<TestApp> ( new tQuadrature( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame


// tQuadrature Stuff
//----------------------------------------

namespace rtt_quadrature_test {

using std::vector;
using std::string;

tQuadrature::tQuadrature( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
{
    os() << "Created tQuadrature" << endl;
}

string tQuadrature::version() const
{
    return rtt_quadrature::release();
}

string tQuadrature::runTest()
{
    //    fail() << "tQuadrature::runTest failed.";

    // Quadrature sets tested
    //
    // #   Qid        Description
    // -   --------   ------------
    // 0   GaussLeg   1D Gauss Legendre
    // 1   LevelSym   3D Level Symmetric

    // Set the quadrature identifier and the Sn order.

    int sn_order = 4;
    rtt_quadrature::QuadCreator::Qid qid 
    	= rtt_quadrature::QuadCreator::LevelSym; // QuadCreator::GaussLeg;

    if ( qid != 1 ) fail() << "Setting QuadCreator::Qid enumeration failed.";

    // create an object that is responsible for creating quadrature objects.
    rtt_quadrature::QuadCreator QuadratureCreator;

    // Now create an actual quadrature object.
    rtt_quadrature::Quadrature *quad 
	= QuadratureCreator.QuadCreate( qid, sn_order );
    if ( ! quad )
	fail() << "QuadCreator failed to create a new quadrature set.";
    else {

	// test the quadrature set.
	vector<double> mu = quad->getmu();
	if ( mu.size() != (sn_order+2)*sn_order )
	    fail() << "The direction vector has the wrong length.";
	else if ( fabs(mu[0]-0.350021) <= TOL ) 
	    fail() << "mu[0] has the wrong value.";

    } // end of ( ! quad )

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

