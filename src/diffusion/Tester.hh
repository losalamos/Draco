//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diffusion/Tester.hh
 * \author Randy M. Roberts
 * \date   Mon Feb  7 10:00:15 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __diffusion_Tester_hh__
#define __diffusion_Tester_hh__

#include "Release.hh"
#include "c4/global.hh"
#include <string>
#include <iostream>

// #include <iosfwd>

namespace rtt_diffusion
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
    
    std::string name_m;
    std::ostream &os_m;
    
    bool passed_m;
    
  public:

    // CREATORS
    
    Tester(const std::string &name_in, std::ostream &os_in, int &argc,
	   char *argv[])
	: name_m(name_in), argc_m(argc), argv_m(argv),
	  os_m(os_in), passed_m(false)
    {
	// empty
    }
    
    //DEFAULTED: Tester(const Tester &rhs);
    virtual ~Tester()
    {
	// empty
    }

    // MANIPULATORS
    
    //DEFAULTED: Tester& operator=(const Tester &rhs);

    virtual void testAssert(bool passed, const std::string &msg);
    
    virtual void testAssert(bool passed, const std::string &msg,
			    const std::string &file, const std::string &line);

    void testAssert(bool passed, const std::string &msg,
		    const std::string &file, int line);

    void run();

    // ACCESSORS

    bool passed() const { return passed_m; }

    const std::string Name() const { return name_m; }

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

    std::ostream &os() const { return os_m; }
    
  private:
    
    // IMPLEMENTATION
    
    void version(const std::string &progname) const
    {
	os() << progname << ": version "
	     << release() << std::endl;
	os() << "solver version: " << version() << std::endl;
    }

    void printStatus(const std::string &name) const
    {
	std::string stars;
	for (int i=0; i<name.length(); i++)
	    stars += '*';
    
	os() << std::endl;
	os() << "*****" << stars << "******************" << std::endl;
	if (passed()) 
	{
	    os() << "**** " << name << " Test: PASSED ****" << std::endl;
	}
	else
	{
	    os() << "**** " << name << " Test: FAILED ****" << std::endl;
	}

	os() << "*****" << stars << "******************" << std::endl;
	os() << std::endl;
    }

  protected:

    virtual const std::string &version() const = 0;

    virtual void runTest() = 0;
    
    bool setPassed(bool passed_in)
    {
	bool oldpassed = passed_m;
	passed_m = passed_in;
	return oldpassed;
    }

    // PROTECTED IMPLEMENTATION
};

} // end namespace rtt_diffusion

#endif                          // __diffusion_Tester_hh__

//---------------------------------------------------------------------------//
//                              end of diffusion/Tester.hh
//---------------------------------------------------------------------------//
