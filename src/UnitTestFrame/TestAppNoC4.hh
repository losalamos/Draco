//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/TestAppNoC4.hh
 * \author Randy M. Roberts
 * \date   Fri Feb 25 07:34:16 2000
 * \brief  The abstract base class, derived from TestAppBase, for non-C4 usage.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __UnitTestFrame_TestAppNoC4_hh__
#define __UnitTestFrame_TestAppNoC4_hh__

#include "TestAppBase.hh"

/*!
 * The only magic here is that ScalarProcessPolicy.hh is included,
 * so that non-C4 versions of ProcessInit and ProcessFinalize
 * are defined in the rtt_UnitTestFrame namespace.
 */
 
#include "ScalarProcessPolicy.hh"

namespace rtt_UnitTestFrame
{

/*!
 * \class TestAppNoC4
 * \brief Abstract base class for test applications, for use with C4.
 */

struct TestAppNoC4 : public TestAppBase
{
    static rtt_dsxx::SP<TestAppNoC4> create(int &argc, char *argv[],
                                            std::ostream &os_in);
  public:
    
    TestAppNoC4(int &argc,  char *argv[], std::ostream &os_in = std::cerr)
        : TestAppBase(argc, argv, os_in)
    {
        // Empty
    }
};

rtt_dsxx::SP<TestAppBase> TestAppBase::create(int &argc, char *argv[],
                                              std::ostream &os_in)
{
    return rtt_dsxx::SP<TestAppBase>(TestAppNoC4::create(argc, argv, os_in));
}

// Define a typedef to TestApp so test programs can always use the
// same name, i.e. TestApp.

typedef TestAppNoC4 TestApp;

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_TestAppNoC4_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/TestAppNoC4.hh
//---------------------------------------------------------------------------//
