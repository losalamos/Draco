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

#include <list>
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
class UnStructuredMeshCFAVTest 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef typename MTFactory::MT MT;
    typedef typename MTFactory::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

  public:

    //! MsgList typedef is a list of pass/fail flags
    //! and their associated messages.

    typedef std::list< std::pair<bool, std::string> > MsgList;

  private:
    
    // DATA
    
    MTFactory &meshFactory_m;

    MsgList msgList_m;
    
  public:

    // CREATORS

    //! Constructor
    
    UnStructuredMeshCFAVTest(MTFactory &meshFactory_in)
	: meshFactory_m(meshFactory_in)
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
	UnStructuredMeshCFAVTest::TestsNotYetWritten();
    }

    // ACCESSORS

    //! Returns a list of pass/fail flags, and their associated messages.

    const MsgList &msgList() const { return msgList_m; }

  private:
    
    // DISSALLOWED CREATORS

    UnStructuredMeshCFAVTest(const UnStructuredMeshCFAVTest &rhs);

    // DISSALLOWED MANIPULATORS
    
    UnStructuredMeshCFAVTest& operator=(const UnStructuredMeshCFAVTest &rhs);

    // IMPLEMENTATION

    void addMsg(bool passed, const std::string &msg)
    {
	msgList_m.push_back(std::pair<bool, std::string>(passed, msg));
    }
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_UnStructuredMeshCFAVTest_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/UnStructuredMeshCFAVTest.hh
//---------------------------------------------------------------------------//
