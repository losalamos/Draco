//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/PassFailStream.hh
 * \author Randy M. Roberts
 * \date   Fri Feb 25 07:58:58 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __UnitTestFrame_PassFailStream_hh__
#define __UnitTestFrame_PassFailStream_hh__

#include "TestApp.hh"
#include <iosfwd>
#include <sstream>

namespace rtt_UnitTestFrame
{
 
//===========================================================================//
/*!
 * \class PassFailStreamBase
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Ch, class Tr>
class PassFailStreamBase
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    std::ostringstream ost_m;
    
  public:

    // CREATORS
    
    PassFailStreamBase()
    {
	// empty
    }
    
    //Defaulted: PassFailStreamBase(const PassFailStreamBase &rhs);

    //Defaulted: ~PassFailStreamBase()

    // MANIPULATORS
    
    //Defaulted: PassFailStreamBase& operator=(const PassFailStreamBase &rhs);

    template<class T>
    PassFailStreamBase &operator<<(const T& val)
    {
	ost() << val;
	return *this;
    }

    // Handle Manipulators
    
    PassFailStreamBase &operator<<(std::basic_ostream<Ch,Tr>
				   &(*f)(std::basic_ostream<Ch,Tr> &))
    {
	f(ost());
	return *this;
    }

    PassFailStreamBase &operator<<(std::ios_base &(*f)(std::ios_base &))
    {
	f(ost());
	return *this;
    }

    PassFailStreamBase &operator<<(std::basic_ios<Ch,Tr>
				   &(*f)(std::basic_ios<Ch,Tr> &))
    {
	f(ost());
	return *this;
    }

    // ACCESSORS

  protected:

    // PROTECTED MANIPULATORS

    std::ostringstream &ost() { return ost_m; }
   
  private:
    
    // IMPLEMENTATION
};

class PassFailStream : public PassFailStreamBase<char, std::char_traits<char> >
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    TestApp &testApp;
    const bool stateToSet;

  public:

    // CREATORS
    
    PassFailStream(TestApp &testApp_in, const std::string &str,
		   bool stateToSet_in)
	: PassFailStreamBase<char, std::char_traits<char> >(),
	  testApp(testApp_in), stateToSet(stateToSet_in)
    {
	ost() << "**** " << testApp.name();
	if (str.size() != 0)
	{
	    ost() << "(" << str << ")";
	}
	ost() << " Test: "
	    << (stateToSet ? "PASSED" : "FAILED")
	    << ": ****\n\t";

	if (stateToSet == false)
	    testApp.setPassed(false);
    }

    //Defaulted: PassFailStream(const PassFailStream &rhs);

    ~PassFailStream()
    {
	testApp.addMessage(ost().str());
    }
    
    // MANIPULATORS
    
    //Defaulted: PassFailStream& operator=(const PassFailStream &rhs);

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_PassFailStream_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/PassFailStream.hh
//---------------------------------------------------------------------------//
