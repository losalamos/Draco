//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Assert.hh
 * \author Geoffrey Furnish, Kelly Thompson
 * \date   6 December 1993
 * \brief  Header file for Draco specific exception class definition
 *         (rtt_dsxx::assertion). Also define Design-by-Contract macros.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_dsxx_Assert_HH
#define RTT_dsxx_Assert_HH

#include <stdexcept>
#include <string>

// add ds++ package configure
#include <ds++/config.h>

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class assertion Exception notification class for Draco specific
 * assertions.
 *
 * This class is derived from std::runtime_error.  In fact, this class
 * provides no significant change in functionality from std::runtime_error.
 * This class provides the following features in addition to those found in
 * std::runtime_error: 
 *
 * 1. rtt_dsxx::assertion does provide an alternate constructor that allows
 *    us to automatically insert file name and line location into error
 *    messages.  
 *
 * 2. It provides a specialized form of std::runtime_error.  This allows
 *    Draco code to handle Draco specific assertions differently from generic
 *    C++ or STL exceptions.  For example
 *
 * \code
 *    try
 *    {
 *       throw rtt_dsxx::assertion( "My error message." );
 *    } 
 *    catch ( rtt_dsxx::assertion &a ) 
 *    {
 *       // Catch Draco exceptions first.
 *       cout << a.what() << endl;
 *       exit(1);
 *    }
 *    catch ( std::runtime_error &e )
 *    { 
 *       // Catch general runtime_errors next
 *       cout << e.what() << endl;
 *    }
 *    catch ( ... )
 *    {
 *       // Catch everything else
 *        exit(1);
 *    }
 * \endcode
 *
 * \note Assertion should always be thrown as objects on the stack and caught
 *       as references. 
 */
/*!
 * \example ds++/test/tstAssert.cc
 * 
 * Assertion and DBC examples.
 */
// revision history:
// -----------------
// 11 March 2003 (Kelly Thompson)
// 
// Major changes for class rtt_dsxx::assertion:
//
// Let rtt_dsxx::assertion derive from std::runtime_error.  G. Furnish
// intended for the design to follow this model but the C++ compilers at time
// did not have full support for <stdexcept>.  API is unchanged but the guts
// of the class are significantly different (i.e.: use copy and assignment
// operators from base class, use string instead of const char *, etc.).
//
// 28 July 1997
//
// Total rewrite of the DS++ Assertion facility.  In the old formulation,
// Assert and Insist expanded into inline if statements, for which the body
// were blocks which both contained output to cerr, and which threw an
// exception.  This has proven problematic because it meant that Assert.hh
// included iostream.h, and also because throw was inlined (not a problem for
// Assert per se, since Assert was generally shut off for optimized codes,
// but still a problem for Insist).  These effects combined to present
// problems both at compile time and, possibly, at run time.
//
// The new formulation, motivated by suggestions from Arch Robison and Dave
// Nystrom, uses a call to an out-of-line function to do the actual logging
// and throwing of exceptions.  This permits removing iostream.h from the
// inclusion graph for translation units including Assert.hh, and avoids
// inlining of exception throwing code.  Moreover, it also allows us to get
// __FILE__ into play, which we couldn't do before because there was no way
// to do string concatenation with it, and you couldn't do free store
// management effectively in an inlined macro.  By going out to a global
// function, we pick up the ability to formulate a more complete picture of
// the error, and provide some optimization capability (both compile time and
// run time).
//
// Note also that at this juncture, we go ahead and drop all support for
// compilers which are incapable of compiling exception code.  From here
// forward, exceptions are assumed to be available.  And thereby we do our
// part to promote "standard C++".
//===========================================================================//

class assertion : public std::logic_error
{
  public:
    /*!
     * \brief Default constructor for ds++/assertion class.
     *
     * This constructor creates a ds++ exception object.  This object is
     * derived form std::runtime_error and has identical functionality.  The
     * principal purpose for this class is to provide an exception class that
     * is specialized for Draco.  See the notes on the overall class for more
     * details.
     *
     * \param msg The error message saved with the exception.
     */
    explicit assertion( std::string const & msg )
	:  std::logic_error( msg )
    { /* empty */ }

    /*!
     * \brief Specialized constructor for rtt_dsxx::assertion class.
     *
     * This constructor creates a ds++ exception object.  This object is
     * derived form std::runtime_error and has identical functionality.  This
     * constructor is specialized for use by Draco DbC commands (Require,
     * Ensure, Check, and Insist).  It forms the error message from the test
     * condition and the file and line number of the DbC command.
     *
     * \param cond The expression that failed a DbC test.
     * \param file The source code file name that contains the DbC test.
     * \param line The source code line number that contains the DbC test.
     *
     * \sa \ref Draco_DBC, --with-dbc[=level], Require, Ensure, Check, Insist
     */
    assertion( std::string const & cond, 
	       std::string const & file, 
	       int const line )
	: std::logic_error( build_message( cond, file, line ) )
    { /* empty */ }

    /*!
     * \brief Destructor for ds++/assertion class.
     *
     * We do not allow the destructor to throw!
     */
    virtual ~assertion() throw() { /* empty */ }

  private:
    /*!
     * \brief Helper function to build error message that includes source
     * file name and line number.
     */
    std::string build_message( std::string const & cond, 
			       std::string const & file, 
			       int         const   line ) const;
};

//---------------------------------------------------------------------------//
// FREE NAMESPACE FUNCTIONS
//---------------------------------------------------------------------------//

// Throw a rtt_dsxx::assertion for Require, Check, Ensure.
void toss_cookies( std::string const & cond, 
		   std::string const & file, 
		   int         const line );

// Throw a rtt_dsxx::assertion for Insist.
void insist( std::string const & cond, 
	     std::string const & msg, 
	     std::string const & file, 
	     int         const line);

} // end of rtt_dsxx

//---------------------------------------------------------------------------//
/*!
 * \page Draco_DBC Using the Draco Design-by-Contract Macros
 *
 * The assertion macros are intended to be used for validating preconditions
 * which must be true in order for following code to be correct, etc.  For
 * example, Assert( x > 0. ); y = sqrt(x); If the assertion fails, the code
 * should just bomb.  Philosophically, it should be used to feret out bugs
 * in preceding code, making sure that prior results are within reasonable
 * bounds before proceeding to use those results in further computation,
 * etc.
 * 
 * These macros are provided to support the Design By Contract formalism.
 * The activation of each macro is keyed off a bit in the DBC macro which can 
 * be specified on the command line:
 *    Bit     DBC macro affected
 *     0       Require
 *     1       Check
 *     2       Ensure
 * So for instance, -DDBC=7 turns them all on, -DDBC=0 turns them all
 * off, and -DDBC=1 turns on Require but turns off Check and Ensure.  The
 * default is to have them all enabled.
 *
 * The Insist macro is akin to the Assert macro, but it provides the
 * opportunity to specify an instructive message.  The idea here is that you
 * should use Insist for checking things which are more or less under user
 * control.  If the user makes a poor choice, we "insist" that it be
 * corrected, providing a corrective hint.
 * 
 * \note We provide a way to eliminate assertions, but not insistings.  The
 * idea is that Assert is used to perform sanity checks during program
 * development, which you might want to eliminate during production runs for
 * performance sake.  Insist is used for things which really really must be
 * true, such as "the file must've been opened", etc.  So, use Assert for
 * things which you want taken out of production codes (like, the check might
 * inhibit inlining or something like that), but use Insist for those things
 * you want checked even in a production code.
 */
/*!
 * \def Require(condition)
 * 
 * Pre-condition checking macro.  On when DBC & 1 is true.
 */
/*!
 * \def Check(condition)
 * 
 * Intra-scope checking macro.  On when DBC & 2 is true.
 */
/*!
 * \def Ensure(condition)
 * 
 * Post-condition checking macro.  On when DBC & 4 is true.
 */
/*!
 * \def Remember(code)
 * 
 * Add code to compilable code.  Used in the following manner:
 * \code
 *     Remember (int old = x;)
 *     // ...
 *     Ensure (x == old);
 * \endcode
 * On when DBC & 4 is true.
 */
/*!
 * \def Insist(condition, message)
 * 
 * Inviolate check macro.  Insist is always on.
 */
//---------------------------------------------------------------------------//

#if !defined(DBC)
#define DBC 7
#endif

#if DBC & 1
#define Require(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ );
#else
#define Require(c) 
#endif

#if DBC & 2
#define Check(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ );
#define Assert(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ );
#else
#define Check(c) 
#define Assert(c) 
#endif

#if DBC & 4
#define Ensure(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ );
#define Remember(c) c;
#else
#define Ensure(c) 
#define Remember(c)
#endif

#define Insist(c,m) if (!(c)) rtt_dsxx::insist( #c, m, __FILE__, __LINE__ );

#endif				// RTT_dsxx_Assert_HH

//---------------------------------------------------------------------------//
//                              end of ds++/Assert.hh
//---------------------------------------------------------------------------//
