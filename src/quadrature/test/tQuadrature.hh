//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Quadrature/tQuadrature.hh
 * \author Kelly Thompson
 * \date   Mon Mar 6 13:41:03 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __Quadrature_test_tQuadrature_hh__
#define __Quadrature_test_tQuadrature_hh__

#include "UnitTestFrame/TestApp.hh"

namespace rtt_quadrature_test
{
 
//===========================================================================//
/*!
 * \class tQuadrature
 *
 * \brief A class used to test the QuadCreator and Quadrature classes.  It is 
 * inheireted from Randy's UnitTestFrame::TestApp.
 */
// revision history:
// -----------------
// 0) original
// 1) 3-17-2000 a) Added comments (some in Doxygen format).
//              b) Changed the test program to loop over a sequence of
//                 quadrature types performing the same tests on each
//                 object. 
//              c) Implemented use of ds++/Assert class.
// 
//===========================================================================//

class tQuadrature : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    tQuadrature( int argc, char *argv[], std::ostream &os_in );
    //Defaulted:  tQuadrature(const tQuadrature &rhs);
    //Defaulted: ~tQuadrature();

    // MANIPULATORS
    
    //Defaulted: tQuadrature& operator=(const tQuadrature &rhs);

    // ACCESSORS

    std::string name() const { return "tQuadrature"; }
    std::string version() const;

  protected:
    
    std::string runTest();

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_Quadrature_test

#endif // __Quadrature_test_tQuadrature_hh__

//---------------------------------------------------------------------------//
//		      end of Quadrature/tQuadrature.hh
//---------------------------------------------------------------------------//
