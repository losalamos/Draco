//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfException.hh
 * \author Kelly Thompson
 * \date   Tue Sep  5 10:47:29 2000
 * \brief  GandolfException class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfException_hh__
#define __cdi_gandolf_GandolfException_hh__

#include <string>
#include <stdexcept>

namespace rtt_cdi_gandolf
{
 
    //========================================================================
    /*!
     * \class GandolfException
     *
     * \brief This class handles exceptions thrown by GandolfOpacity when
     *        it calls the Gandolf library functions.
     *
     * This class provides an Gandolf exception data object derived
     * from std::exception.  Additionally, there are Gandolf exception
     * objects (derived from GandolfException) for each function that
     * Gandolf provides.  There is one base exception class
     * (GandolfException) and 5 derived exception classes
     * (gkeysException, gchgridsException, ggetgrayException,
     * ggetmgException and gmatidsException).
     *
     * When an exception is thrown in the cdi_gandolf package a
     * XXXException (see list above for what XXX might stand for)
     * object is created and may be trapped in the calling routine by
     * using a try/catch block when a GandolfFile or GandolfOpacity
     * object is first created or when the client fetches an opacity
     * value(s) by using one of GandolfOpacities accessor functions.
     *
     * \example cdi_gandolf/test/tGandolfFile.cc
     * \example cdi_gandolf/test/tCDIGandolf.cc
     * \example cdi_gandolf/test/tGandolfWithCDI.cc
     *
     * These examples demonstrate how the try/catch blocks may be used
     * to trap errors in the cdi_gandolf package.  They also
     * demonstrate how the calling program can extract information
     * about the exception.
     *  
     */
    //========================================================================

    class GandolfException : public std::exception
    { 
      private:

	/*!
	 * \brief The name of the Gandolf library function that
	 *	 returned an error.
	 */
	const std::string gandolfFunctionName;

	/*!
	 * \brief A summary message that describes the error that was
	 *	 caught.
	 */
	mutable char* message;
	
      protected:
	
	/*!
	 * \brief The Gandolf error code that was returned from the
	 *	 FORTRAN library.
	 */
 	const int errorCode;
	
      public:

	// CREATORS

	/*!
	 * \brief The standard GandolfExceptoin constuctor.
	 *
	 * When a specific error is thrown (i.e.: gmatidException),
	 * the base class is also instantiated through this 
	 * constructor.  This constructor, in turn, instantiates 
	 * std::exception.
	 *
	 * The throw() at the end of these declarations are exception
	 * specifiers.  I have not limited the types of exceptions
	 * that may be caught by this object (hense the empty parens).
	 * See Lippman and Lajoie, "C++ Primer" p. 1042ff for more
	 * information.
  	 *
	 * \param what_arg A simple description of the exception.  In
	 *        some cases this may simply be the name of the
	 *        Gandolf function that failed.  In other cases it
	 *        might contain a detailed description of the error.
	 *        This string is used in the construction of
	 *        std::exception.
	 * \param gandolfFunctionName This optional argument is a
	 *        string holds the name of the Gandolf function that
	 *        failed.
	 * \param errorCode This optional argument is the integer
	 *        error code returned by Gandolf.  GandolfExcetion can
	 *        translate this into a string message. 
	 */
	GandolfException( const std::string& what_arg,
			  std::string gandolfFunctionName = "unknown",
			  int errorCode = -99 ) throw();
	
	/*!
	 * \brief The standard GandolfException destructor.
	 *
	 * \sa what()
	 *
	 * The destructor frees memory that may have been allocated by
	 * the accessor what().  If what was never called then
	 * "message" was never allocated and no deallocation is
	 * required.
	 */
	virtual ~GandolfException() throw();
	
	// ACCESSORS

	/*!
	 * \brief GandolfException overrides the default
	 *        std::exception.what() function so that a more
	 *        detailed message may be returned. 
	 * 
	 * \sa The actual text message returned is generated by
	 *     getErrorSummary(). 
	 * \sa ~GandolfException(): Memory allocated by this routine
	 *     is freed in the class destrcutor.
	 *
	 * This function returns the contents of the data member
	 * "message".  If message is empty a call to errorSummary() is
	 * made to fill "message" with descriptive text.  If required
	 * memory is allocated for the data member "message".
	 * Deallocation is completed in the destructor
	 * ~GandolfException().
	 */
	virtual const char* what() const throw();

	/*!
	 * \brief Returns the numeric error code supplied by the
	 *        Gandolf function. 
	 */
	int getErrorCode() const { return errorCode; };

	/*!
	 * \brief Returns a brief message that is correlated to the
	 * 	  specific error code returned by Gandolf.
	 *
	 * This is a virtual function that is overridden by the
	 * derived Gandolf exception classes.  Each of derived classes
	 * can decode the numeric Gandolf error code into a short
	 * descriptive string.
	 */
	virtual std::string getErrorMessage() const;

	/*!
	 * \brief Returns the name of the Gandolf function that
	 *	  returned an error.
	 */
	const std::string& getGandolfFunctionName() const 
	{ 
	    return gandolfFunctionName; 
	};

	/*!
	 * \brief Returns a detailed message containing the Gandolf
	 *        function that failed, the actual error code returned
	 *        and a brief description of the error.
	 */
	std::string errorSummary() const;
    };


    //========================================================================
    /*!
     * \class gkeysException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function gkeys().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     *
     * There are two possible constructors.  One simply takes a user
     * supplied error message.  The second takes the integer supplied
     * Gandolf error code and looks up the associated error message.
     *
     * getErrorMessage() overrides the base class and provides the
     * ability to lookup an integer error code and find its associated
     * error message.
     */
    //========================================================================

    class gkeysException : public GandolfException
    {
      public:
	
	// CREATORS

	gkeysException( const std::string& what_arg );
	gkeysException( int gkeysError );
	~gkeysException() throw();
	virtual std::string getErrorMessage() const;
    };
       
   //========================================================================
    /*!
     * \class gchgridsException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function gchgrids().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     *
     * There are two possible constructors.  One simply takes a user
     * supplied error message.  The second takes the integer supplied
     * Gandolf error code and looks up the associated error message.
     *
     * getErrorMessage() overrides the base class and provides the
     * ability to lookup an integer error code and find its associated
     * error message.
     */
    //========================================================================

    class gchgridsException : public GandolfException
    {
      public:
	gchgridsException( const std::string& what_arg );
	gchgridsException( int gchgridsError );
	~gchgridsException() throw();
	virtual std::string getErrorMessage() const;
    };

    //========================================================================
    /*!
     * \class ggetgrayException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function ggetgray().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     *
     * There are two possible constructors.  One simply takes a user
     * supplied error message.  The second takes the integer supplied
     * Gandolf error code and looks up the associated error message.
     *
     * getErrorMessage() overrides the base class and provides the
     * ability to lookup an integer error code and find its associated
     * error message.
     */
    //========================================================================

    class ggetgrayException : public GandolfException
    {
      public:
	ggetgrayException( const std::string& what_arg );
	ggetgrayException( int ggetgrayError );
	~ggetgrayException() throw ();
	virtual std::string getErrorMessage() const;
    };

    //========================================================================
    /*!
     * \class ggetmgException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function ggetmg().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     *
     * There are two possible constructors.  One simply takes a user
     * supplied error message.  The second takes the integer supplied
     * Gandolf error code and looks up the associated error message.
     *
     * getErrorMessage() overrides the base class and provides the
     * ability to lookup an integer error code and find its associated
     * error message.
     */
    //========================================================================

     class ggetmgException : public GandolfException
    {
      public:
	ggetmgException( const std::string& what_arg );
	ggetmgException( int ggetmgError );
	~ggetmgException() throw();
	virtual std::string getErrorMessage() const;
    };

    //========================================================================
    /*!
     * \class gmatidsException
     *
     * \brief Derived from GandolfException this class encapsulates
     *        errors returned by the Gandolf Function gmatids().
     * 
     * \sa The class description for GandolfException for additional
     *     comments.
     *
     * There are two possible constructors.  One simply takes a user
     * supplied error message.  The second takes the integer supplied
     * Gandolf error code and looks up the associated error message.
     *
     * getErrorMessage() overrides the base class and provides the
     * ability to lookup an integer error code and find its associated
     * error message.
     */
    //========================================================================

    class gmatidsException : public GandolfException
    {
      public:
	gmatidsException( const std::string& what_arg );
	gmatidsException( int gmatidsError );
	~gmatidsException() throw();
	virtual std::string getErrorMessage() const;
    };

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfException_hh__

//---------------------------------------------------------------------------//
// end of cdi_gandolf/GandolfException.hh
//---------------------------------------------------------------------------//
