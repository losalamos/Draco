//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/TestApp.hh
 * \author Randy M. Roberts
 * \date   Fri Feb 25 07:34:16 2000
 * \brief  The abstract base class, derived from TestAppBase, for C4 usage.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __UnitTestFrame_TestApp_hh__
#define __UnitTestFrame_TestApp_hh__

#include "TestAppBase.hh"

/*!
 * The only magic here is that C4ProcessPolicy.hh is included,
 * so that C4 versions of ProcessInit and ProcessFinalize
 * are defined in the rtt_UnitTestFrame namespace.
 */
 
#include "C4ProcessPolicy.hh"

namespace rtt_UnitTestFrame
{

/*!
 * \class TestApp
 * \brief Abstract base class for test applications, for use with C4.
 */

struct TestApp : public TestAppBase
{
    static rtt_dsxx::SP<TestApp> create(int &argc, char *argv[],
                                        std::ostream &os_in);
  public:
    
    TestApp(int &argc,  char *argv[], std::ostream &os_in = std::cerr)
        : TestAppBase(argc, argv, os_in)
    {
        // Empty
    }
};

rtt_dsxx::SP<TestAppBase> TestAppBase::create(int &argc, char *argv[],
                                              std::ostream &os_in)
{
    return rtt_dsxx::SP<TestAppBase>(TestApp::create(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_TestApp_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/TestApp.hh
//---------------------------------------------------------------------------//
