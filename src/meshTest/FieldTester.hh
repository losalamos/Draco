//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/FieldTester.hh
 * \author Shawn Pautz, Randy M. Roberts
 * \date   Mon Aug 23 16:49:40 1999
 * \brief  Header file for the FieldTester class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_FieldTester_hh__
#define __meshTest_FieldTester_hh__

#include "Tester.hh"

#include <iostream>
#include <string>

namespace rtt_meshTest
{
 
//===========================================================================//
/*!
 * \class FieldTester
 * \brief Class used by TestMTFields in order to test compliance
 *  of MT nested containers to the concepts needed by "Solon."
 */
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, class FGRP>
class FieldTester : public Tester
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MT::FieldConstructor FieldConstructor;

    typedef typename FGRP::XD XD;
    typedef typename FGRP::XI XI;
    typedef typename FGRP::XL XL;
    typedef typename FGRP::XDC XDC;

    // DATA

    FieldConstructor &fCtor_m;
    MT &mesh_m;
    
  public:

    // CREATORS
    
    FieldTester(const std::string &name_in,
		FieldConstructor &fCtor_in, MT &mesh_in,
		std::ostream &os_in)
	: Tester(std::string("FieldTester<") + name_in + ">", os_in),
	  fCtor_m(fCtor_in), mesh_m(mesh_in)
    {
	/* empty */
    }
	
    ~FieldTester() { /* empty */ };

    // MANIPULATORS
    
    void run();

    // ACCESSORS

  private:
    
    // DISSALLOWED CREATORS

    FieldTester(const FieldTester &rhs);

    // DISSALLOWED MANIPULATORS
    
    FieldTester& operator=(const FieldTester &rhs);

    // IMPLEMENTATION


    void t1();
    void t2();
    void t3();
    void t4();
    void t5();
    void t6();
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_FieldTester_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/FieldTester.hh
//---------------------------------------------------------------------------//
