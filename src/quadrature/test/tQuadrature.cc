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

namespace rtt_quadrature_test {

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::fabs;
using rtt_dsxx::SP;
using rtt_quadrature::QuadCreator;
using rtt_quadrature::Quadrature;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
/*!
 * \brief Tests the Quadcrator and Quadtrature constructors and access
 * routines. 
 *
 * To add a quadrature to this test the following items must be changed:
 *   add new enumeration to Qid[] array.
 *   add new mu[0] value to mu0[] array.
 *   verify nquads is set to the correct number of quadrature sets being
 *      tested. 
 */
void quadrature_test()
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
    
    SP< const Quadrature > spQuad;

    // loop over quadrature types to be tested.

    for ( int ix = 0; ix < nquads; ++ix ) {
	
	// Verify that the enumeration value matches its int value.
	if ( qid[ix] != ix ) {
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
	    cout << "   Number of Angles = " << spQuad->getNumAngles() << endl;

	    // If the object was constructed sucessfully then we continue
	    // with the tests.
	    if ( ! spQuad )
		FAILMSG("QuadCreator failed to create a new quadrature set.")
	    else {
		// get the mu vector
		vector<double> mu = spQuad->getMu();
		// get the omega vector for direction m=1.
		vector<double> omega_1 = spQuad->getOmega(1);
		// get the total omega object.
		vector< vector<double> > omega = spQuad->getOmega();
		// compare values.
		if ( mu.size() != spQuad->getNumAngles() )
		    FAILMSG("The direction vector has the wrong length.")
		else if ( fabs( mu[0] + mu0[ix] ) >= TOL ) 
		    FAILMSG("mu[0] has the wrong value.")
		else if ( fabs( mu[1] - omega_1[0] ) >= TOL )
		    FAILMSG("mu[1] != omega_1[0].")
		else if ( fabs( omega[1][0] - omega_1[0] ) >= TOL )
		    FAILMSG("omega[1][0] != omega_1[0].")
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
    return;
} // end of quadrature_test


//===========================================================================//
// PASS/FAILURE
//===========================================================================//

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

bool fail(int line, char *file)
{
    std::cout << "Test: failed on line " << line << " in " << file
	      << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

bool pass_msg(const std::string &passmsg)
{
    std::cout << "Test: passed" << std::endl;
    std::cout << "     " << passmsg << std::endl;
    return true;
}

//---------------------------------------------------------------------------//

bool fail_msg(const std::string &failmsg)
{
    std::cout << "Test: failed" << std::endl;
    std::cout << "     " << failmsg << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//
// BOOLEAN PASS FLAG
//---------------------------------------------------------------------------//

bool passed = true;

} // end namespace rtt_quadrature_test


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;
    using std::string;

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_quadrature::release()
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	rtt_quadrature_test::quadrature_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tQuadrature, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_quadrature_test::passed) 
    {
        cout << "**** tQuadrature Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tQuadrature." << endl;

    return 0;
}   



//---------------------------------------------------------------------------//
//                            end of tQuadrature.cc
//---------------------------------------------------------------------------//

