//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/tQuadServices.cc
 * \author Kelly Thompson
 * \date   Mon Nov 8 10:48 2004
 * \brief  Quadrature Services.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
//#include <cmath>
//#include <sstream>
#include <string>

#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include "quadrature_test.hh"
#include "../Quadrature.hh"
#include "../QuadCreator.hh"
#include "../QuadServices.hh"
#include "../Release.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// //---------------------------------------------------------------------------//
// /*! 
//  * \brief test_legendre_poly
//  * 
//  * Test a function that evaluates legendre polynomials.
//  */
// void test_legendre_poly()
// {
//     using rtt_dsxx::soft_equiv;
//     using rtt_dsxx::SP;
//     using rtt_quadrature::QuadCreator;
//     using rtt_quadrature::Quadrature;
//     using rtt_quadrature::QuadServices;

//     // Create a quadrature set from a temporary instance of a
//     // QuadratureCreator factory object.
//     SP< const Quadrature > spQuad;
//     spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym, 2 ); 

//     // Create a QuadServices object so we can test its member functions.
//     QuadServices qs( spQuad );

//     // P_{0,0}(x) == 1.0
//     double reference( 1.0 );
//     unsigned k(0), ell(0);
//     double value( qs.legendre_polynomial( k, ell, 1.0 ) );

//     if( soft_equiv(reference,value) )
//     {
// 	PASSMSG("Correct computation of P00(x).");
//     }
//     else
//     {
// 	FAILMSG("Failed to compute P00(x).");
//     }
    
//     return;
// }


//---------------------------------------------------------------------------//
/*!
 * \brief Tests the Quadcrator and Quadtrature constructors and access
 * routines. 
 *
 * To add a quadrature to this test the following items must be changed: add
 * new enumeration to Qid[] array.  add new mu[0] value to mu0[] array.
 * verify nquads is set to the correct number of quadrature sets being
 * tested.
 */
void quadrature_test()
{
     using rtt_quadrature::QuadCreator;
     using rtt_quadrature::Quadrature;
     using rtt_quadrature::QuadServices;
     using rtt_dsxx::SP;
 //     using rtt_dsxx::soft_equiv;
    
     using std::cout;
     using std::endl;
     using std::string;
     using std::vector;

     //----------------------------------------
     // Setup Quadrature set

     // create an object that is responsible for creating quadrature objects.
     // QuadCreator QuadratureCreator;
    
     // we will only look at S2 Sets in this test.
     const size_t sn_ord_ref( 2                    );
     const string qname_ref ( "3D Level Symmetric" );
     const size_t n_ang_ref ( 8                    );

     // Banner
     cout << "Testing the "  << qname_ref << " quadrature set." << endl;

     // Create a quadrature set from a temporary instance of a
     // QuadratureCreator factory object.
     SP< const Quadrature > spQuad;
     spQuad = QuadCreator().quadCreate( QuadCreator::LevelSym, sn_ord_ref ); 
    
     // print the name of the quadrature set that we are testing.
     const string qname   ( spQuad->name()         );
     const size_t sn_order( spQuad->getSnOrder()   );
     const size_t n_ang   ( spQuad->getNumAngles() );

     // check basic quadrature setup.
     if( sn_order != sn_ord_ref ) ITFAILS;
     if( n_ang    != n_ang_ref  ) ITFAILS;
     if( qname    != qname_ref  ) ITFAILS;

     // Print a table
     spQuad->display();

     //----------------------------------------
     // Setup QuadServices object

     QuadServices qs( spQuad );

     vector<double> const M( qs.getM() );
     unsigned const numMoments( M.size()/n_ang );
     cout << "iang momt value" << endl
	  << "---- ---- --------------------------------" <<endl;
     for( unsigned m=0; m<n_ang; ++m )
	 for( unsigned n=0; n<numMoments; ++n )
	     cout << "  " << m << "    " << n << "   "
		  << M[ n + m*numMoments ] << endl;


    return;
} // end of quadrature_test

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::string;
    using std::cout;
    using std::endl;

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
// 	test_legendre_poly();
	quadrature_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tQuadServices, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_quadrature_test::passed) 
    {
        cout << "**** tQuadServices Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tQuadServices." << endl;

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tQuadServices.cc
//---------------------------------------------------------------------------//
