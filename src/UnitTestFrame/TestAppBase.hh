//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/TestAppBase.hh
 * \author Randy M. Roberts
 * \date   Fri Feb 25 07:34:16 2000
 * \brief  The abstract base from which unit test applications can be derived.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __UnitTestFrame_TestAppBase_hh__
#define __UnitTestFrame_TestAppBase_hh__

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
 * \class TestAppBase
 * This class forms the abstract base from which unit test applications
 * can be derived.  Unit tests derived from this class must supply the
 * name() and version() methods, as well as the runTest() method, which will
 * be called to perform the actual unit test.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestAppBase 
{

    // NESTED CLASSES AND TYPEDEFS AND FRIEND CLASSES

    /*!
     * The PassFailStream class makes use of several protected methods
     * of TestAppBase objects, inorder to perform its functions.
     */
    
    friend class PassFailStream;

    // DATA

    /*!
     * The static theSPTestAppBase smart pointer is an attempt to enforce
     * the Singleton pattern.
     */
    
    static rtt_dsxx::SP<TestAppBase> theSPTestAppBase;

    /*!
     * The argument list, defined by argc and argv in the constructor call,
     * are stored in a list of std::string's.
     */
    
    std::list<std::string> argList_m;

    /*!
     * This message list contains all of the pass/fail messages.
     * Pass and fail messages are not sent directly to the ostream,
     * but are stored in the message list for later processing.
     */
    
    std::list<std::string> messageList_m;

    /*!
     * A flag to determine whether this unit test passed or failed.
     */
    
    bool passed_m;

    /*!
     * The ostream to which pass/fail messages should be written.
     */
    
    std::ostream &os_m;

  public:

    // CREATORS

    /*!
     * The ctor for the unit test base class.  The argc and argv represent
     * the command line arguments.  They will be used to create a list of
     * command line argument std::string's.  The os_in argument refers to
     * the ostream to which pass/fail messages are written.
     */
    
    TestAppBase(int &argc,  char *argv[], std::ostream &os_in = std::cerr);
    
    //Defaulted: TestAppBase(const TestAppBase &rhs);

    /*!
     * The dtor is virtual, since derived class' dtor's must be called
     * via base class references or pointers.
     */
    
    virtual ~TestAppBase()
    {
        // Empty
    }

    // MANIPULATORS
    
    //Defaulted: TestAppBase& operator=(const TestAppBase &rhs);

    /*!
     * The run method is the public interface used to run the unit test(s).
     * The run method, in turn, calls the protected virtual runTest() method.
     * If you are using the mainDriver's main(argc,argv) then this method
     * is called automatically (outside of any try blocks).
     * The return value is printed to the pass/fail message ostream.
     * If --version is supplied in the arguments, run() returns the version
     * info and does not run runTest.
     */
    
    std::string run();

    // PROTECTED MANIPULATORS
    
  protected:

    /*!
     * The run method is the protected virtual method used to run the unit
     * test(s).
     * This method is called within the public run() method.
     * The derived class must supply its own version of this class.
     * The return value is usually returned by run().
     */
    
    virtual std::string runTest() = 0;

    /*!
     * Get your grubby little hands on the argument list.
     */
    
    std::list<std::string> &argList() { return argList_m; }

    /*!
     *             !!!YOUR MESSAGE HERE!!!
     *
     * Add your own messages, just like the pass/fail messages
     * added by the PassFailStream objects.  If you are using the
     * mainDriver framework then these messages will be printed
     * to the ostream.
     */
    
    void addMessage(const std::string &message)
    {
	messageList_m.push_back(message);
    }

    /*!
     * Protected method to change the pass/fail status of the test.
     */
    
    bool setPassed(bool flag)
    {
	bool tmp = passed_m;
	passed_m = flag;
	return tmp;
    }
    
    /*!
     * Now the magic...
     * The pass() method is for use within the overloaded runTest() methods.
     * This method returns a PassFailStream just aching to have data output
     * onto them via the operator<< "inserter" (Stroustrup page 608).
     * The pass() method gives you a PassFailStream that says a subtest passed.
     */
    
    PassFailStream pass(const std::string &str = "");

    /*!
     * The fail() method is for use within the overloaded runTest() methods.
     * This method returns a PassFailStream just aching to have data output
     * onto them via the operator<< "inserter" (Stroustrup page 608).
     * The fail() method gives you a PassFailStream that says a subtest failed.
     */
    
    PassFailStream fail(const std::string &str = "");

  public:

    // ACCESSORS

    /*!
     * The ostream to which pass/fail messages should be written.
     */
    
    std::ostream &os() const { return os_m; }

    /*!
     * The pass/fail state of this unit test.
     */
    
    bool passed() const { return passed_m; }

    /*!
     * The pass/fail message list.
     */
    
    const std::list<std::string> &messageList() const { return messageList_m; }

    /*!
     * A convenience utility that retrieves the argument after the specified
     * argument, i.e. foo --a x --b y --z --t
     * with a call to getNextArg("--b") returns "y".
     */
    
    std::string getNextArg(const std::string &arg) const;

    /*!
     * The name of the unit test.  Must be supplied by derived class.
     */
    
    virtual std::string name() const = 0;

    /*!
     * The version of the unit being tested.  Must be supplied by the
     * derived class.
     */
    
    virtual std::string version() const = 0;

  public:
    
    // STATIC CREATOR

    /*!
     * The public static method of accessing the **one** object of the class
     * derived from this class.  Used to enforce the Singleton pattern.
     * This method calls the static private create method.
     */
    
    static TestAppBase &theTestAppBase(int &argc, char *argv[],
			       std::ostream &os_in = std::cerr)
    {
	if (!theSPTestAppBase)
        {
           theSPTestAppBase = create(argc, argv, os_in);
        }
	return *theSPTestAppBase;
    }

    /*!
     * The public static method of accessing the **one** object of the class
     * derived from this class.  Used to enforce the Singleton pattern.
     * This method calls the static private create method.
     */
    
    static TestAppBase &theTestAppBase()
    {
	int argc = 0;
	char **argv = 0;
	if (!theSPTestAppBase)
	    theSPTestAppBase = create(argc, argv, std::cerr);
	return *theSPTestAppBase;
    }

    // STATIC CLASS METHODS

    /*!
     * The public static method that initializes the process,
     * creates a singleton TestAppBase derived object,
     * and calls the run() method on it.
     */

    static int driveTest(int &argc, char *argv[]);

  private:

    // STATIC IMPLEMENTATION
    
    /*!
     * When setting up a unit test the author of the test must **define**
     * this static method that creates a new object of the derived class,
     * and returns a base class smart pointer to it.
     * For Example, if MyTestApp is derived from TestAppBase,
     * then in MyTestApp.cc (ignoring namespaces)...

     * SP<TestAppBase> TestAppBase::create(int &argc, char *argv[],
     *                                     ostream& os_in)
     * {
     *    return SP<TestAppBase>(new MyTestApp(argc, argv, os_in));
     * }
     */
    
     static rtt_dsxx::SP<TestAppBase> create(int &argc, char *argv[],
                                             std::ostream &os_in);
    
  private:

    // IMPLEMENTATION

};

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_TestAppBase_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/TestAppBase.hh
//---------------------------------------------------------------------------//
