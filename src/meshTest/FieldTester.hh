//----------------------------------*-C++-*----------------------------------//
// FieldTester.hh
// Randy M. Roberts
// Mon Aug 23 16:49:40 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __meshTest_FieldTester_hh__
#define __meshTest_FieldTester_hh__

#include <string>

namespace rtt_meshTest
{
 
//===========================================================================//
// class FieldTester - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, class FGRP>
class FieldTester 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MT::FieldConstructor FieldConstructor;

    typedef typename FGRP::XD XD;
    typedef typename FGRP::XI XI;
    typedef typename FGRP::XL XL;
    typedef typename FGRP::XDC XDC;

    // DATA

    std::string name_m;
    
    FieldConstructor &fCtor_m;
    MT &mesh_m;
    
    std::ostream &os_m;
    
    bool passed_m;

  public:

    // CREATORS
    
    FieldTester(const std::string &name_in,
		FieldConstructor &fCtor_in, MT &mesh_in,
		std::ostream &os_in)
	: name_m(name_in), fCtor_m(fCtor_in), mesh_m(mesh_in),
	  os_m(os_in), passed_m(false)
    {
	/* empty */
    }
	
    ~FieldTester() { /* empty */ };

    // MANIPULATORS
    
    void run();

    // ACCESSORS

    bool passed() const { return passed_m; }

  private:
    
    // DISSALLOWED CREATORS

    FieldTester(const FieldTester &rhs);

    // DISSALLOWED MANIPULATORS
    
    FieldTester& operator=(const FieldTester &rhs);

    // IMPLEMENTATION

    void error(bool &passed, const std::string &msg);

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
