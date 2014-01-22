//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/ApplicationUnitTest.hh
 * \author Kelly Thompson
 * \date   Thu Jun  1 17:15:05 2006
 * \brief  Declaration file for encapsulation of Draco unit test for applications. 
 * \note   Copyright (C) 2006-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
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

#include "ds++/UnitTest.hh"
#include "Timer.hh"

#if defined(MSVC)
#   pragma warning (push)
#   pragma warning (disable:4251) //  warning C4251: 'rtt_c4::ApplicationUnitTest::applicationName' : class 'std::basic_string<_Elem,_Traits,_Ax>' needs to have dll-interface to be used by clients of class 'rtt_c4::ApplicationUnitTest'
#endif

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
 * \warning The output from this class is closely tied to the CTest
 * Pass/Fail regular expressions listed in \c config/component_macros.cmake
 * and in the local \c CMakeLists.txt. 
 *
 */
/*! 
 * \example c4/test/tstApplicationUnitTest.cc 
 * This unit test demonstrates typical usage for ApplicationUnitTest. * 
 */
//===========================================================================//

class DLL_PUBLIC ApplicationUnitTest : public rtt_dsxx::UnitTest
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
    // Each successful return code increments numPasses.  Each error return
    // code increments numFails.  So a successful series of tests
    // (listOfArgs.size()>1) will have at least this many passes and no
    // fails.  numPasses can be greater than listOfArgs.size() if the  user
    // has incremented the pass count manually.
    bool allTestsPass(void) const
    { return numPasses>=listOfArgs.size() && numFails == 0; };
    
    //! Provide a report of the number of unit test passes and fails.
    void status(void);

    //! Name of the log file
    std::string logFileName() const { return logFile; }

    //! Return number of processors that are to be used by this test.
    std::string nodes() const { return numProcs; }

    //! Return whether problem computation times should be reported
    bool reportTimingsI() const { return reportTimings; }
        
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

    //! Should the times for each test be reported?
    bool reportTimings;

    //! Timer if timings are desired
    Timer problemTimer;
};

} // end namespace rtt_c4

#if defined(MSVC)
#   pragma warning (pop)
#endif

#endif // c4_ApplicationUnitTest_hh

//---------------------------------------------------------------------------//
// end of c4/ApplicationUnitTest.hh
//---------------------------------------------------------------------------//
