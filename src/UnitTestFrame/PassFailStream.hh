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
#include <sstream>

namespace rtt_UnitTestFrame
{
 
//===========================================================================//
/*!
 * \class PassFailStream
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class PassFailStream
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    TestApp &testApp;
    const bool stateToSet;
    std::ostringstream ost;
    
  public:

    // CREATORS
    
    PassFailStream(TestApp &testApp_in, const std::string &str,
		   bool stateToSet_in)
	: testApp(testApp_in), stateToSet(stateToSet_in)
    {
	ost << "**** " << testApp.name();
	if (str.size() != 0)
	{
	    ost << "(" << str << ")";
	}
	ost << " Test: "
		<< (stateToSet ? "PASSED" : "FAILED")
		<< ": ****\n\t";

	if (stateToSet == false)
	    testApp.setPassed(false);
    }
    
    ~PassFailStream()
    {
	testApp.addMessage(ost.str());
    }

    // MANIPULATORS
    
    template<class T>
    PassFailStream &operator<<(const T& val)
    {
	ost << val;
	return *this;
    }

    // ACCESSORS

  private:
    
    // DISALLOWED METHODS
    
    PassFailStream(const PassFailStream &rhs);
    PassFailStream& operator=(const PassFailStream &rhs);

    // IMPLEMENTATION
};

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_PassFailStream_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/PassFailStream.hh
//---------------------------------------------------------------------------//
