//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/UnStructuredMeshCFAVTest.hh
 * \author Randy M. Roberts
 * \date   Wed Sep 22 13:23:07 1999
 * \brief  Header file for the base class UnStructuredMeshCFAVTest
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_UnStructuredMeshCFAVTest_hh__
#define __meshTest_UnStructuredMeshCFAVTest_hh__

#include "Tester.hh"
#include <string>

namespace rtt_meshTest
{
 
//===========================================================================//
/*!
 * \class UnStructuredMeshCFAVTest
 * \brief Class used by TestMTConnFacesArroundVrtx to test
 * the unstructured MT::ConnFacesAroundVertices service.
 *
 * The UnStructuredMeshCFAVTest class is used by TestMTConnFacesArroundVrtx
 * to run tests for an
 * unstructured Mesh MT.  These tests excersize the
 * MT::ConnFacesAroundVertices service.
 *
 */
// 
//===========================================================================//

template<class MTFactory>
class UnStructuredMeshCFAVTest  : public Tester
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

  private:
    
    // DATA

    Tester &parent_m;
    
    MTFactory &meshFactory_m;

  public:

    // CREATORS

    //! Constructor
    
    UnStructuredMeshCFAVTest(Tester &parent_in, MTFactory &meshFactory_in)
	: parent_m(parent_in), meshFactory_m(meshFactory_in)
    {
	/* empty */
    }

    //! Destructor
    
    ~UnStructuredMeshCFAVTest() { /* empty */ }

    // MANIPULATORS

    //! The run method is the primary interface into the class.
    //! It has not yet been implemented, and is designed to fail
    //! at compile time.
    
    void run()
    {
	setPassed(true);
	UnStructuredMeshCFAVTest::TestsNotYetWritten();
    }

    // ACCESSORS

    // PROTECTED MANIPULATORS

    virtual void testassert(bool passed, const std::string &msg)
    {
	parent_m.testassert(passed, Name()+": "+msg);
	if (!passed)
	    setPassed(false);
    }
    
    virtual void testassert(bool passed, const std::string &msg,
			    const std::string &file,
			    const std::string &line)
    {
	Tester::testassert(passed, msg, file, line);
    }
    
    virtual void testassert(bool passed, const std::string &msg,
			    const std::string &file, int line)
    {
	Tester::testassert(passed, msg, file,line);
    }
    
  private:

    // DISSALLOWED CREATORS

    UnStructuredMeshCFAVTest(const UnStructuredMeshCFAVTest &rhs);

    // DISSALLOWED MANIPULATORS
    
    UnStructuredMeshCFAVTest& operator=(const UnStructuredMeshCFAVTest &rhs);

    // IMPLEMENTATION
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_UnStructuredMeshCFAVTest_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/UnStructuredMeshCFAVTest.hh
//---------------------------------------------------------------------------//
