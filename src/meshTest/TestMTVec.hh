//----------------------------------*-C++-*----------------------------------//
// TestMTVec.hh
// Randy M. Roberts
// Fri Aug 20 09:11:42 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __meshTest_TestMTVec_hh__
#define __meshTest_TestMTVec_hh__

#include <iosfwd>
#include <string>

namespace rtt_meshTest
{
 
//===========================================================================//
// class TestMTVec - 
//
// Purpose :  Test the MT::ccvsf::value_type container according to the Random
//            Access Container requirements.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MTFactory>
class TestMTVec 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

    typedef typename MT::ccvsf ccvsf;
    typedef typename ccvsf::value_type XVEC;

   // DATA

    MTFactory &meshFactory_m;
    
    std::ostream &os_m;
    
    bool passed_m;

  public:

    // CREATORS
    
    TestMTVec(MTFactory &meshFactory_in, std::ostream &os_in)
	: meshFactory_m(meshFactory_in), os_m(os_in),
	  passed_m(false)
    {
	/* empty */
    }
    
    ~TestMTVec() { /* empty */ }

    // MANIPULATORS

    void run();

    // ACCESSORS

    bool passed() const { return passed_m; }
    
  private:
    
    // DISSALLOWED CREATORS

    TestMTVec(const TestMTVec &rhs);

    // DISSALLOWED MANIPULATORS
    
    TestMTVec& operator=(const TestMTVec &rhs);

    // IMPLEMENTATION

    void error(bool &passed, const std::string &msg);

    void t1();
    void t2();
    void t3();
    void t4();
    void t5();
    void t6();
    void t7();
    void t8();
    
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_TestMTVec_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/TestMTVec.hh
//---------------------------------------------------------------------------//
