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

#include "TestAppBase.hh"
#include <iosfwd>
#include <sstream>

namespace rtt_UnitTestFrame
{
 
//===========================================================================//
/*!
 * \class PassFailStreamBase
 *
 * The base class for PassFailStream's.
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

    /*!
     * operator to stream out a data to the pass/fail stream
     */
    
    template<class T>
    PassFailStreamBase &operator<<(const T& val)
    {
	ost() << val;
	return *this;
    }

    // Handle Manipulators
    
    /*!
     * operator to stream out a iostream manipulators to the pass/fail stream
     */
    
    PassFailStreamBase &operator<<(std::basic_ostream<Ch,Tr>
				   &(*f)(std::basic_ostream<Ch,Tr> &))
    {
	f(ost());
	return *this;
    }

    /*!
     * operator to stream out a iostream manipulators to the pass/fail stream
     */
    
    PassFailStreamBase &operator<<(std::ios_base &(*f)(std::ios_base &))
    {
	f(ost());
	return *this;
    }

    /*!
     * operator to stream out a iostream manipulators to the pass/fail stream
     */
    
    PassFailStreamBase &operator<<(std::basic_ios<Ch,Tr>
				   &(*f)(std::basic_ios<Ch,Tr> &))
    {
	f(ost());
	return *this;
    }

    // ACCESSORS

  protected:

    // PROTECTED MANIPULATORS

    /*!
     * Protected access to the ostringstream representation.
     * There must be some way to not expose the representation to
     * derived classes, but I haven't devoted enough energy to
     * find one.  Randy M. Roberts
     */
    
    std::ostringstream &ost() { return ost_m; }
   
  private:
    
    // IMPLEMENTATION
};

//===========================================================================//
/*!
 * \class PassFailStream
 *
 * The class that is used to queue up properly formatted pass/fail messages
 * into a TestAppBase object.
 * They messages then can be accessed, presumably to send to std::cerr
 * for parsing.
 */
//===========================================================================//

class PassFailStream : public PassFailStreamBase<char, std::char_traits<char> >
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    TestAppBase &testApp;
    const bool stateToSet;

  public:

    // CREATORS

    /*!
     * This Ctor takes the TestAppBase object reference for which
     * pass/fail messages will be added.
     * The Ctor takes a string that will be added to the pass/fail
     * message in a properly formatted manner.
     * The Ctor takes a boolean value to determine whether this is
     * an accumulator for passing messages, or failing messages.
     * (true -> passing, failse -> failing)
     */
    
    PassFailStream(TestAppBase &testApp_in, const std::string &str,
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

    /*!
     * The dtor adds the accumulated messages to the testApp.
     * With proper usage the pass/fail stream has only the lifetime
     * of one message; therefore, after every message the testApp
     * will be updated.
     */
    
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
