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

#include <iosfwd>
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
    
    std::ostream &os() { return os_m; }
    
    virtual void testassert(bool passed, const std::string &msg);
    
    virtual void testassert(bool passed, const std::string &msg,
			    const std::string &file, const std::string &line);

    void testassert(bool passed, const std::string &msg,
		    const std::string &file, int line);

    // ACCESSORS

    bool passed() const { return passed_m; }

    const std::string Name() const { return name_m; }

    // PROTECTED MANIPULATORS
    
  protected:

    bool setPassed(bool passed_in)
    {
	bool oldpassed = passed_m;
	passed_m = passed_in;
	return oldpassed;
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_Tester_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/Tester.hh
//---------------------------------------------------------------------------//
