//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/TestApp.hh
 * \author Randy M. Roberts
 * \date   Fri Feb 25 07:34:16 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __UnitTestFrame_TestApp_hh__
#define __UnitTestFrame_TestApp_hh__

#include "ds++/SP.hh"

#include <string>
#include <list>
#include <iostream>

namespace rtt_UnitTestFrame
{

// Forward Reference

class PassFailStream;

//===========================================================================//
/*!
 * \class TestApp
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestApp 
{

    // NESTED CLASSES AND TYPEDEFS AND FRIEND CLASSES

    friend class PassFailStream;

    // DATA

    static rtt_dsxx::SP<TestApp> theSPTestApp;

    std::list<std::string> argList_m;
    std::list<std::string> messageList_m;
    
    bool passed_m;
    std::ostream &os_m;

  public:

    // CREATORS
    
    TestApp(int &argc,  char *argv[], std::ostream &os_in = std::cerr);
    
    //Defaulted: TestApp(const TestApp &rhs);

    virtual ~TestApp() { /* empty */ }

    // MANIPULATORS
    
    //Defaulted: TestApp& operator=(const TestApp &rhs);

    std::string run();

    // PROTECTED MANIPULATORS
    
  protected:

    virtual std::string runTest() = 0;

    std::list<std::string> &argList() { return argList_m; }
    void addMessage(const std::string &message)
    {
	messageList_m.push_back(message);
    }
    
    PassFailStream pass(const std::string &str = "");
    PassFailStream fail(const std::string &str = "");

  public:

    // ACCESSORS

    std::ostream &os() const { return os_m; }

    bool passed() const { return passed_m; }

    const std::list<std::string> &messageList() const { return messageList_m; }

    virtual std::string name() const = 0;

    virtual std::string version() const = 0;

  public:
    
    // STATIC CREATOR

    static TestApp &theTestApp(int &argc, char *argv[])
    {
	if (!theSPTestApp)
	    theSPTestApp = create(argc, argv);
	return *theSPTestApp;
    }

    static TestApp &theTestApp()
    {
	int argc = 0;
	char **argv = 0;
	if (!theSPTestApp)
	    theSPTestApp = create(argc, argv);
	return *theSPTestApp;
    }

  private:

    // STATIC IMPLEMENTATION
    
    static rtt_dsxx::SP<TestApp> create(int &argc, char *argv[]);

  private:

    // IMPLEMENTATION

    bool setPassed(bool flag)
    {
	bool tmp = passed_m;
	passed_m = flag;
	return tmp;
    }
    
};

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_TestApp_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/TestApp.hh
//---------------------------------------------------------------------------//
