//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/Tester.hh
 * \author Randy M. Roberts
 * \date   Mon Nov 22 16:08:52 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_Tester_hh__
#define __meshTest_Tester_hh__

#include <iostream>
#include <sstream>
#include <string>

namespace rtt_meshTest
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

    std::string name_m;
    
    std::ostream &os_m;
    
    bool passed_m;
    
  public:

    // CREATORS

    Tester(const std::string &name_in,
	   std::ostream &os_in)
	: name_m(name_in), os_m(os_in), passed_m(false)
    {
	// empty
    }
    
    virtual ~Tester() { };

    // MANIPULATORS
    
    void testassert(bool passed, const std::string &msg)
    {
	if (!passed)
	{
	    os_m << name_m << " failed: " << msg << std::endl;
	    passed_m = false;
	}
    }


    void testassert(bool passed, const std::string &msg,
		    const std::string &file, const std::string &line)
    {
	testassert(passed, msg + " in " + file + " line: " + line);
    }

    void testassert(bool passed, const std::string &msg,
		    const std::string &file, int line)
    {
	std::ostringstream strline;
	strline << line;
	testassert(passed, msg, file, strline.str());
    }

    // ACCESSORS

    bool passed() const { return passed_m; }

    // PROTECTED MANIPULATORS
    
  protected:

    bool setPassed(bool passed_in)
    {
	bool oldpassed = passed_m;
	passed_m = passed_in;
	return oldpassed;
    }

    std::ostream &os() { return os_m; }
    
  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_Tester_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/Tester.hh
//---------------------------------------------------------------------------//
