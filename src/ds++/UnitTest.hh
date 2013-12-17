//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/UnitTest.hh
 * \author Kelly Thompson
 * \date   Thu May 18 15:46:19 2006
 * \brief  Provide some common functions for unit testing within Draco
 * \note   Copyright (C) 2006-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_UnitTest_hh
#define dsxx_UnitTest_hh

#include "Assert.hh"
#include <list>
#include <iostream>
#include <map>
#include <vector>

#if defined(MSVC)
#   pragma warning (push)
#   pragma warning (disable:4251) //  warning C4251: 'rtt_dsxx::UnitTest::testName' : class 'std::basic_string<_Elem,_Traits,_Ax>' needs to have dll-interface to be used by clients of class 'rtt_dsxx::UnitTest'
#endif

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class UnitTest
 * \brief Object to encapsulate unit testing of Draco classes and functions. 
 *
 * This is a virtual class.  You should use one of the following UnitTest
 * classes in your test appliation:
 *
 * \li ScalarUnitTest      - Used for testing code that does not use parallel
 *                       communication (rtt_c4).
 * \li ParallelUnitTest    - Used for testing code that does use parallel
 *                       communications (rtt_c4). 
 * \li ApplicationUnitTest - Used for testing applications that run in
 *                       parallel. The unit test code is run in scalar-mode
 *                       and calls mpirun to run the specified application.
 *
 *
 * \sa UnitTest.cc for additional details.
 *
 * \par Code Sample:
 *
 * Scalar UnitTests should have the following syntax.
 * \code

#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)
#define ITFAILS    ut.failure( __LINE__, __FILE__ )
 
int main(int argc, char *argv[])
{
    try
    {
        // Test ctor for ScalarUnitTest (also tests UnitTest ctor)
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
 * Or you can use the UT_EPILOG to shorten the test main() function:
 *
 * \code
int main(int argc, char *argv[])
{
    rtt_utils::ScalarUnitTest ut( argc, argv, release );
    try
    {
        // Test ctor for ScalarUnitTest (also tests UnitTest ctor)
        tstOne(ut);
    }
    UT_EPILOG(ut);
}
 * \endcode
 *
 * \test All of the member functions of this class are tested by
 * ds++/test/tstScalarUnitTest.cc, including the early exit caused by 
 * \c --version on the command line.
 */
//===========================================================================//

class DLL_PUBLIC UnitTest 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    //! Typedef for function pointer to this package's release function.
    typedef std::string const (*string_fp_void)(void);
    
    // CREATORS
        
    //! Default constructors.
    UnitTest( int    & argc,
              char **& argv,
              string_fp_void   release_,
              std::ostream   & out_         = std::cout );
    
    //! The copy constructor is disabled.
    UnitTest( UnitTest const &rhs );

    //! Destructor is virtual because this class will be inherited from.
    virtual ~UnitTest(void){/*empty*/};

    // MANIPULATORS    

    //! The assignment operator is disabled.
    UnitTest& operator=( UnitTest const &rhs );

    //! Only special cases should use these (like the unit test
    //! tstScalarUnitTest.cc). 
    void dbcRequire( bool b ) { m_dbcRequire=b; return; }
    void dbcCheck(   bool b ) { m_dbcCheck=b;   return; }
    void dbcEnsure(  bool b ) { m_dbcEnsure=b;  return; }

    //! Change the target for output
    // void setostream( std::ostream out_ ) { out = out_; return; };
    
    // ACCESSORS
    bool failure(int line);
    bool failure(int line, char const *file);
    bool failure( std::string const &failmsg );
    bool passes(  std::string const &passmsg );
    //! This pure virtual function must be provided by the inherited class.
    //It should provide output concerning the status of UnitTest.
    void status(void) const { out << resultMessage() << std::endl; return; }
    //! Reset the pass and fail counts to zero.
    void reset() { numPasses=0;numFails=0; return; }
    bool dbcRequire(void) const { return m_dbcRequire; }
    bool dbcCheck(void)   const { return m_dbcCheck;   }
    bool dbcEnsure(void)  const { return m_dbcEnsure;  }
    bool dbcNothrow(void) const { return m_dbcNothrow; }
    bool dbcOn(void)      const { return m_dbcRequire || m_dbcCheck || m_dbcEnsure; }
    std::string getTestPath(void) const { return testPath; }
    std::string getTestName(void) const { return testName; }
    std::string getTestInputPath(void) const;

    // DATA
    //! The number of passes found for this test.
    unsigned numPasses;
    //! The number of failures found for this test.
    unsigned numFails;

    //! Is fpe_traping active?
    bool fpe_trap_active;

    // Features
    static std::map< std::string, unsigned >
    get_word_count( std::ostringstream const & data, bool verbose=false );
    static std::map< std::string, unsigned >
    get_word_count( std::string const & filename, bool verbose=false );
    static std::vector<std::string> tokenize(
        std::string const & source,
        char        const * delimiter_list = " ",
        bool                keepEmpty      = false);
    
  protected:

    // IMPLEMENTATION
    std::string resultMessage(void) const;

    // DATA
    
    //! The name of this unit test.
    std::string const testName;
    //! Relative path to the unit test.
    std::string const testPath;

    //! Function pointer this package's release(void) function
    string_fp_void release;

    //! Where should output be sent (default is std::cout)
    std::ostream & out;

    /*! Save the state of DBC so that it is easily accessible from within a
     * unit test.
     */
    bool m_dbcRequire;
    bool m_dbcCheck;
    bool m_dbcEnsure;
    bool m_dbcNothrow;

};

} // end namespace rtt_dsxx

//#define PASSMSG(ut,m)  ut.passes(m)
//#define FAILMSG(ut,m)  ut.failure(m)
//#define ITFAILS(ut)    ut.failure( __LINE__, __FILE__ )
//#define FAILURE(ut)    ut.failure( __LINE__, __FILE__ );
//#define UT_PROLOG(foo) typedef ut foo
#define UT_EPILOG(foo) \
catch (rtt_dsxx::assertion &err) {     \
   std::cout << "DRACO ERROR: While testing " << foo.getTestName() << ", " \
             << "the following error was thrown...\n" \
             << err.what() << std::endl; foo.numFails++; } \
catch(std::exception &err) { \
   std::cout << "ERROR: While testing " << foo.getTestName() << ", " \
             << "the following error was thrown...\n" \
             << err.what() << std::endl; foo.numFails++; } \
catch( ... ) { \
   std::cout << "ERROR: While testing " << foo.getTestName() << ", " \
             << "An unknown exception was thrown on processor " \
             << std::endl; foo.numFails++; }; \
return foo.numFails; 

#if defined(MSVC)
#   pragma warning (pop)
#endif

#endif // dsxx_UnitTest_hh

//---------------------------------------------------------------------------//
// end of ds++/UnitTest.hh
//---------------------------------------------------------------------------//
