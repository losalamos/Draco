//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/ParallelUnitTest.hh
 * \author Kelly Thompson
 * \date   Thu Jun  1 17:15:05 2006
 * \brief  Declaration file for encapsulation of Draco parallel unit tests.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 *
 * This file provides a definition for ParallelUnitTest.  The purpose of this
 * class is to encapsulate the keywords and behavior of DBS parallel unit
 * tests. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_ParallelUnitTest_hh
#define c4_ParallelUnitTest_hh

#include "SpinLock.hh"
#include "global.hh"
#include "ds++/UnitTest.hh"
#include <iostream>

namespace rtt_c4
{

//===========================================================================//
/*!
 * \class ParallelUnitTest
 * \brief This class encapsulates services for parallel unit tests.
 *
 * This class inherits from UnitTest.  Much of the documentation for the
 * services of this class is provided in UnitTest.hh
 *
 * \sa rtt_dsxx::UnitTest for detailed description of all the UnitTest
 * classes.
 *
 * \par Code Sample:
 * \code
int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        tstOne(ut);
        tstTwo(ut);
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstSwap, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstSwap, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
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
 * changes to these DBS m4 macro command
 */
/*! 
 * \example c4/test/tstParallelUnitTest.cc 
 * This unit test demonstrates typical usage for ParallelUnitTest. * 
 */
//===========================================================================//

class ParallelUnitTest : public rtt_dsxx::UnitTest
{
  public:

    // CREATORS
    
    //! Default constructor.
    ParallelUnitTest( int            & argc,
                      char         **& argv,
                      string_fp_void   release_,
                      std::ostream   & out_ = std::cout);

    //!  The copy constructor is disabled.
    ParallelUnitTest( ParallelUnitTest const &rhs );

    //! Destructor.
    ~ParallelUnitTest();

    // MANIPULATORS
    
    //! The assignment operator is disabled.
    ParallelUnitTest& operator=( ParallelUnitTest const &rhs );

    // ACCESSORS
    
    //! Provide a report of the number of unit test passes and fails.
    void status(void);
    
  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA

};

} // end namespace rtt_c4

#endif // c4_ParallelUnitTest_hh

//---------------------------------------------------------------------------//
//              end of c4/ParallelUnitTest.hh
//---------------------------------------------------------------------------//
