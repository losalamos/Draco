//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/test/Tester.hh
 * \author Randy M. Roberts
 * \date   Wed Jan 26 13:01:54 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMG_test_Tester_hh__
#define __LAMG_test_Tester_hh__

#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"

#include <string>
#include <iostream>

#define TESTASSERT(c) testAssert(c, #c, __FILE__, __LINE__)

namespace rtt_LAMG_test
{
 
//===========================================================================//
/*!
 * \class Tester
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Tester 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
    int &argc_m;
    char **argv_m;
    std::string filename_m;
    bool passed_m;
    std::ostream &os_m;
    
  public:

    // CREATORS
    
    Tester(const std::string &filename, int &argc, char **argv,
	   std::ostream &os_in)
	: argc_m(argc), argv_m(argv), filename_m(filename),
	  passed_m(false), os_m(os_in)
    {
	std::cout << "Calling C4::Init(...)" << std::endl;
	
	C4::Init( argc_m, argv_m );

    }
	
    //DEFAULTED: Tester(const Tester &rhs);

    virtual ~Tester()
    {
	std::cout << "Calling C4::Finalize()" << std::endl;
	C4::Finalize();
    }

    // MANIPULATORS
    
    //DEFAULTED: Tester& operator=(const Tester &rhs);

    inline void run();

    // ACCESSORS

    virtual const std::string name() const = 0;

    std::string filename() const { return filename_m; }
    bool passed() const { return passed_m; }
    std::ostream &os() { return os_m; }
    const std::string fullName() const { return name(); }

  protected:

    // PROTECTED IMPLEMENTATION

    inline void testAssert(bool passed, const std::string &what,
			   const std::string &file, int line);

    virtual void runTest() = 0;

  private:
    
    // IMPLEMENTATION
    void version(const std::string &progname) const
    {
	std::cout << progname << ": version "
		  << rtt_LAMG::release() << std::endl;
    }

    bool checkForVersion() const
    {
	for (int arg=1; arg < argc_m; arg++)
	{
	    if (std::string(argv_m[arg]) == "--version")
	    {
		version(argv_m[0]);
		return true;
	    }
	}
	return false;
    }

    void printStatus(const std::string &str) const
    {
	C4::HTSyncSpinLock sl;
	
	std::string stars;
	for (int i=0; i<str.length(); i++)
	    stars += '*';
    
	std::cout << std::endl;
	std::cout << "*****" << stars << "******************" << std::endl;
	if (passed()) 
	{
	    std::cout << "**** " << str << " Test: PASSED ****" << std::endl;
	}
	else
	{
	    std::cout << "**** " << str << " Test: FAILED ****" << std::endl;
	}

	std::cout << "*****" << stars << "******************" << std::endl;
	std::cout << std::endl;
	std::cout << std::flush;
    }
};

void Tester::run()
{
    if (checkForVersion())
	return;
    
    try
    {
	std::cout << "Initiating "
		  << fullName()
		  << " test"
		  << std::endl;

	passed_m = true;
    
	runTest();

    }
    catch( dsxx::assertion& a )
    {
	std::cout << "Failed assertion: " << a.what() << std::endl;
	passed_m = false;
    }
    catch( std::exception& a )
    {
	std::cout << "Failed std::eception: " << a.what() << std::endl;
	passed_m = false;
    }

    // Print the status of the test.

    printStatus(fullName());
}

void Tester::testAssert(bool passed, const std::string &what,
			const std::string &file, int line)
{
    if (!passed)
    {
	os() << fullName()
	     << " failed (" << what << "): in " << file
	     << " line: " << line << "." << std::endl;
	passed_m = false;
    }
}

} // end namespace rtt_LAMG_test

#endif                          // __LAMG_test_Tester_hh__

//---------------------------------------------------------------------------//
//                              end of LAMG/test/Tester.hh
//---------------------------------------------------------------------------//
