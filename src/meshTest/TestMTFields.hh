//----------------------------------*-C++-*----------------------------------//
// TestMTFields.hh
// Randy M. Roberts
// Mon Aug 23 16:49:16 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __meshTest_TestMTFields_hh__
#define __meshTest_TestMTFields_hh__

#include "FieldGroup.hh"
#include "DoubleContainer.hh"

#include <iosfwd>
#include <string>

namespace rtt_meshTest
{
 
//===========================================================================//
// class TestMTFields - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MTFactory>
class TestMTFields 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

  public:

    typedef FieldGroup<typename MT::template bstf<double>,
	typename MT::template bstf<int>, typename MT::template bstf<long>,
	typename MT::template bstf <DoubleContainer> > BSTF;
    
    typedef FieldGroup<typename MT::template cctf<double>,
	typename MT::template cctf<int>, typename MT::template cctf<long>,
	typename MT::template cctf<DoubleContainer> > CCTF;
    
    typedef FieldGroup<typename MT::template fcdtf<double>,
	typename MT::template fcdtf<int>, typename MT::template fcdtf<long>,
	typename MT::template fcdtf<DoubleContainer> > FCDTF;
    
    typedef FieldGroup<typename MT::template nctf<double>,
	typename MT::template nctf<int>, typename MT::template nctf<long>,
	typename MT::template nctf<DoubleContainer> > NCTF;
    
    typedef FieldGroup<typename MT::template vctf<double>,
	typename MT::template vctf<int>, typename MT::template vctf<long>,
	typename MT::template vctf<DoubleContainer> > VCTF;
    
    // DATA

  private:

    MTFactory &meshFactory_m;
    
    std::ostream &os_m;
    
    bool passed_m;

  public:

    // CREATORS
    
    
    TestMTFields(MTFactory &meshFactory_in, std::ostream &os_in)
	: meshFactory_m(meshFactory_in), os_m(os_in),
	  passed_m(false)
    {
	/* empty */
    }
    
    ~TestMTFields() { /* empty */ }

    // MANIPULATORS

    template<class FGRP>
    void run(const std::string &name_in);
    
    // ACCESSORS

    bool passed() const { return passed_m; }

  private:
    
    // DISSALLOWED CREATORS

    TestMTFields(const TestMTFields &rhs);

    // DISSALLOWED MANIPULATORS
    
    TestMTFields& operator=(const TestMTFields &rhs);

    // IMPLEMENTATION
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_TestMTFields_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/TestMTFields.hh
//---------------------------------------------------------------------------//
