//----------------------------------*-C++-*----------------------------------//
// TestMTModel.hh
// Randy M. Roberts
// Fri Aug 20 09:11:42 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __meshTest_TestMTModel_hh__
#define __meshTest_TestMTModel_hh__

#include <iosfwd>

namespace rtt_meshTest
{
 
//===========================================================================//
// class TestMTModel - 
//
// Purpose :  This class tests whether a given MT is a model of the
//            "Solon" MT concept.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MTFactory>
class TestMTModel 
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
    
    TestMTModel(MTFactory &meshFactory_in, std::ostream &os_in)
	: meshFactory_m(meshFactory_in), os_m(os_in),
	  passed_m(false)
    {
	/* empty */
    }
    
    ~TestMTModel() { /* empty */ }

    // MANIPULATORS

    void run();

    // ACCESSORS

    bool passed() const { return passed_m; }
    
  private:
    
    // DISSALLOWED CREATORS

    TestMTModel(const TestMTModel &rhs);

    // DISSALLOWED MANIPULATORS
    
    TestMTModel& operator=(const TestMTModel &rhs);

    // IMPLEMENTATION

    void t1();
    void t2();
    
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_TestMTModel_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/TestMTModel.hh
//---------------------------------------------------------------------------//
