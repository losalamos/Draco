//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/test/Tester.hh
 * \author Randy M. Roberts
 * \date   Mon Nov 22 11:02:57 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_test_Tester_hh__
#define __LAMGDiffusionSolver_test_Tester_hh__

#include "../Release.hh"

#include "utils.hh"

#include "c4/global.hh"

#include <iostream>

#define TESTASSERT(c) testAssert(c, #c, __FILE__, __LINE__)

namespace rtt_LAMGDiffusionSolver_test
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

template<class FACTORY>
class Tester 
{

    // NESTED CLASSES AND TYPEDEFS

  protected:
    
    typedef typename FACTORY::MT MT;
    typedef typename FACTORY::Product MTFactoryProduct;
    typedef typename MT::FieldConstructor FieldConstructor;

    // DATA

  private:

    FACTORY factory_m;
    
    int &argc_m;
    char **argv_m;

    bool passed_m;
    std::ostream &os_m;
    
  public:

    // CREATORS
    
    Tester(const std::string &filename, int &argc, char **argv,
	   std::ostream &os_in)
	: factory_m(getMTFactory(filename)),
	  argc_m(argc), argv_m(argv), passed_m(false), os_m(os_in)
    {
	std::cout << "Calling C4::Init(...)" << std::endl;
	
	C4::Init( argc_m, argv_m );

    }

    virtual ~Tester()
    {
	std::cout << "Calling C4::Finalize()" << std::endl;
	C4::Finalize();
    }

    // MANIPULATORS
    
    inline void run();

    FACTORY &factory() { return factory_m; }

    // ACCESSORS

    const FACTORY &factory() const { return factory_m; }
    bool passed() const { return passed_m; }
    std::ostream &os() { return os_m; }

    virtual const std::string name() const = 0;
    const std::string fullName() const
    {
	return name() + "<" + factory().name() + ">";
    }

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
		  << rtt_LAMGDiffusionSolver::release() << std::endl;
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

template <class FACTORY>
void Tester<FACTORY>::run()
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
    catch( rtt_dsxx::assertion& a )
    {
	std::cout << "Failed assertion: " << a.what() << std::endl;
	passed_m = false;
    }

    // Print the status of the test.

    printStatus(fullName());
}

template<class FACTORY>
void Tester<FACTORY>::testAssert(bool passed, const std::string &what,
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

} // end namespace rtt_LAMGDiffusionSolver_test

#endif                          // __LAMGDiffusionSolver_test_Tester_hh__

//---------------------------------------------------------------------------//
//                              end of LAMGDiffusionSolver/test/Tester.hh
//---------------------------------------------------------------------------//
