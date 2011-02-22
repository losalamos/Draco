//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/ApplicationUnitTest.hh
 * \author Kelly Thompson
 * \date   Thu Jun  1 17:15:05 2006
 * \brief  Declaration file for encapsulation of Draco unit test for applications. 
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 *
 * This file provides a definition for ApplicationUnitTest.  The purpose of
 * this class is to encapsulate the keywords and behavior of DBS application
 * unit tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_ApplicationUnitTest_hh
#define c4_ApplicationUnitTest_hh

#include "c4/config.h"
#include "ds++/UnitTest.hh"
#include <iostream>

namespace rtt_c4
{

//===========================================================================//
/*!
 * \class ApplicationUnitTest
 * \brief This class encapsulates services for application unit tests.
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
    using namespace rtt_c4;
    try
    {
        ApplicationUnitTest ut( argc, argv, release,
                                std::string("../bin/serrano") );
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
 * c4/test/tstApplicationUnitTest.cc, including the early exit caused by 
 * \c --version on the command line.
 *
 * \warning The output from this class is closely tied to the DBS python
 * script \c tools/regression_filter.py that is used during \c gmake \c check.
 * Changing the format or keyword in the output streams from this class should
 * be coordinated with the regular expression matches found in \c
 * tools/regression_filter.py.
 *
 * \warning The differences between ScalarUnitTest, ApplicationUnitTest and
 * ApplicationUnitTest are correlated to the DBS m4 macros \c AC_RUNTESTS and
 * \c AC_TEST_APPLICATION.  Changes to these classes should be coordinated with
 * changes to these DBS m4 macro command
 */
/*! 
 * \example c4/test/tstApplicationUnitTest.cc 
 * This unit test demonstrates typical usage for ApplicationUnitTest. * 
 */
//===========================================================================//

class ApplicationUnitTest : public rtt_dsxx::UnitTest
{
  public:

    // CREATORS
    
    //! Default constructor.
    ApplicationUnitTest(
        int & argc, char **&argv, string_fp_void release_,
        std::string const applicationName_,
        std::list< std::string > const & listOfArgs_=std::list<std::string>(),
        std::ostream & out_ = std::cout ); 

    //! The copy constructor is disabled.
    ApplicationUnitTest(const ApplicationUnitTest &rhs);

    //! Destructor.
    ~ApplicationUnitTest(void){ out << resultMessage() << std::endl; return; };

    // MANIPULATORS

    //! The assignment operator is disabled.
    ApplicationUnitTest& operator=(const ApplicationUnitTest &rhs);

    //! Add a command line argument to the list that will be tested.
    void addCommandLineArgument( std::string const & appArgs );

    //! Execute the specified binary with the provided arguments.
    void runTests(void);
    bool runTest( std::string const & appArgs );

    // ACCESSORS

    //! Did all tests pass?
    bool allTestsPass(void) const
    { return numPasses==listOfArgs.size() && numFails == 0; };
    
    //! Provide a report of the number of unit test passes and fails.
    void status(void);

    //! Name of the log file
    std::string logFileName() const { return logFile; }

    //! Return number of processors that are to be used by this test.
    std::string nodes() const { return numProcs; }
        
  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    //! Construct the MPI directive that will be used to execute the binary.
    std::string constructMpiCommand( std::string const & );

    //! Provide file extension that is used to cpature output from binary.
    std::string buildLogExtension( std::string const & );

    //! Initialize numProcs from command line arguments.
    std::string getNumProcs( int & argc, char **& argv );

    // DATA

    //! The name of the binary file to be tested.
    std::string const applicationName;
    
    //! The path to the binary file to be tested.
    std::string const applicationPath;

    //! Number of processors for the test(s).
    std::string const numProcs;

    //! The command that will be executed by the system.
    std::string const mpiCommand;

    //! The extension used when creating log files
    std::string const logExtension;

    //! A list of command line arguments used during execution of the test.
    std::list< std::string > listOfArgs;

    //! Name of file for logging output of application execution
    std::string logFile;
    
};

} // end namespace rtt_c4


#ifdef c4_isDarwin
#define C4_MPICMD "mpirun -np "
#define C4_UNAME "Darwin"
#endif

#ifdef c4_isLinux
#ifndef C4_MPICMD
#define C4_MPICMD "mpirun -np "
#endif
#define C4_UNAME "Linux"
#endif

#ifdef c4_isOSF1
#define C4_MPICMD "prun -n "
#define C4_UNAME "OSF1"
#endif

#ifdef c4_isAIX
#define C4_MPICMD "poe -procs "
#define C4_UNAME "AIX"
#endif

#ifdef c4_isWin
#define C4_MPICMD "mpiexec -np "
#define C4_UNAME "WINDOWS"
#endif

#ifdef c4_isLinux_with_aprun
#ifdef C4_MPICMD
#undef C4_MPICMD
#endif
#define C4_MPICMD "aprun -n "
#endif

#endif // c4_ApplicationUnitTest_hh

//---------------------------------------------------------------------------//
//              end of c4/ApplicationUnitTest.hh
//---------------------------------------------------------------------------//
