//----------------------------------*-C++-*----------------------------------//
// TestMTComm.hh
// Randy M. Roberts
// Mon Aug 23 15:09:07 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __meshTest_TestMTComm_hh__
#define __meshTest_TestMTComm_hh__

namespace rtt_meshTest
{
 
//===========================================================================//
// class TestMTComm - 
//
// Purpose :This class tests the communications part of a given MT for
//          the "Solon" MT concept.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MTFactory>
class TestMTComm 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

    // DATA
    
    MTFactory &meshFactory_m;
    
    std::ostream &os_m;
    
    bool passed_m;

  public:

    // CREATORS
    
    TestMTComm(MTFactory &meshFactory_in, std::ostream &os_in)
	: meshFactory_m(meshFactory_in), os_m(os_in),
	  passed_m(false)
    {
	/* empty */
    }
    
    ~TestMTComm() { /* empty */ }

    // MANIPULATORS
    
    void run();

    // ACCESSORS

    bool passed() const { return passed_m; }
    
  private:
    
    // DISSALLOWED CREATORS

    TestMTComm(const TestMTComm &rhs);

    // DISSALLOWED MANIPULATORS
    
    TestMTComm& operator=(const TestMTComm &rhs);

    // IMPLEMENTATION

    void t3();
    void t4();
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_TestMTComm_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/TestMTComm.hh
//---------------------------------------------------------------------------//
