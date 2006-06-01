//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/ScalarUnitTest.hh
 * \author Kelly Thompson
 * \date   Thu May 18 17:08:54 2006
 * \brief  Provide services for scalar unit tests
 * \note   Copyright 2006 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_ScalarUnitTest_hh
#define dsxx_ScalarUnitTest_hh

#include <iostream>

#include "UnitTest.hh"

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class ScalarUnitTest
 * \brief This class provides services for scalar unit tests.
 *
 * This class inherits from UnitTest.  Much of the documentation for the
 * services of this class is provided in UnitTest.hh
 *
 * \sa rtt_dsxx::UnitTest for detailed description of all the UnitTest
 * classes.
 *
 * \par Code Sample:
 *
 * Scalar UnitTests should have the following syntax.
 * \code
int main(int argc, char *argv[])
{
    try
    {
        // Test ctor for ScalarUnitTest (also tests UnitTest ctor and member
        // function setTestName).
        rtt_utils::ScalarUnitTest ut( argc, argv, release );
        tstOne(ut);
        ut.status();
    }
    catch( rtt_dsxx::assertion &err )
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        { cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << "An unknown exception was thrown" << endl;
        return 1;
    }
    return 0;
}  
 * \endcode
 *
 * \test All of the member functions of this class are tested by
 * ds++/test/tstScalarUnitTest.cc, including the early exit caused by 
 * \c --version on the command line.
 *
 * \warning The output from this class is closely tied to the DBS python
 * script \c tools/regression_filter.py that is used during \c gmake \c check.
 * Changing the format or keyword in the output streams from this class should
 * be coordinated with the regular expression matches found in \c
 * tools/regression_filter.py.
 *
 * \warning The differences between ScalarUnitTest, ParallelUnitTest and
 * ApplicationUnitTest are correlated to the DBS m4 macros \c AC_RUNTESTS and
 * \c AC_TEST_APPLICATION.  Changes to these classes should be coordinated with
 * changes to these DBS m4 macro commands.
 */
/*!
 * \example ds++/test/tstScalarUnitTest.cc
 * The unit test for and example usage of the ScalarUnitTest class.
 */
//===========================================================================//

class ScalarUnitTest : public UnitTest
{
  public:

    // CREATORS
    
    //! Default constructors.
    ScalarUnitTest( int & argc, char **&argv, string_fp_void release_,
              std::ostream & out_ = std::cout );

    //! The copy constructor is disabled.
    ScalarUnitTest(const ScalarUnitTest &rhs);

    //! Destructor.
    ~ScalarUnitTest(void);

    // MANIPULATORS
    
    //! The assignment operator for ScalarUnitTest is disabled.
    ScalarUnitTest& operator=(const ScalarUnitTest &rhs);

    // ACCESSORS
    void status(void);
    
};

} // end namespace rtt_dsxx

#endif // dsxx_ScalarUnitTest_hh

//---------------------------------------------------------------------------//
//              end of ds++/ScalarUnitTest.hh
//---------------------------------------------------------------------------//
