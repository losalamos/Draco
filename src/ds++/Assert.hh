//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Assert.hh
 * \brief  Header file for Draco specific exception class definition
 *         (rtt_dsxx::assertion). Also define Design-by-Contract macros.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef RTT_dsxx_Assert_HH
#define RTT_dsxx_Assert_HH

#include "ds++/config.h"
#include <stdexcept>
#include <string>

// This warning is issued because we chose to derive from std::logic_error
// (IS-A) instead of creating a class that 'HAS-A' std::logic_error.object.
#if defined(MSVC)
#pragma warning(push)
// non dll-interface class 'std::logic_error' used as base for dll-interface
// class 'rtt_dsxx::assertion'
#pragma warning(disable : 4275)
#endif

namespace rtt_dsxx {

//===========================================================================//
/*!
 * \class assertion
 *
 * Exception notification class for Draco specific assertions.
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
/*!
 * \section assertRevisionHistory Revision history:
 * \subsection date11March2003 11 March 2003 (Kelly Thompson)
 *
 * Major changes for class rtt_dsxx::assertion:
 *
 * Let rtt_dsxx::assertion derive from std::runtime_error.  G. Furnish intended
 * for the design to follow this model but the C++ compilers at time did not
 * have full support for \c <stdexcept>.  API is unchanged but the guts of the
 * class are significantly different (i.e.: use copy and assignment operators
 * from base class, use string instead of \c const \c char \c *, etc.).
 *
 * \subsection date28July1997 28 July 1997
 *
 * Total rewrite of the ds++ Assertion facility.  In the old formulation, Assert
 * and Insist expanded into inline if statements, for which the body were blocks
 * which both contained output to cerr, and which threw an exception.  This has
 * proven problematic because it meant that \c Assert.hh included \c iostream.h,
 * and also because throw was inlined (not a problem for Assert per se, since
 * Assert was generally shut off for optimized codes, but still a problem for
 * Insist).  These effects combined to present problems both at compile time
 * and, possibly, at run time.
 *
 * The new formulation, motivated by suggestions from Arch Robison and Dave
 * Nystrom, uses a call to an out-of-line function to do the actual logging and
 * throwing of exceptions.  This permits removing \c iostream.h from the
 * inclusion graph for translation units including \c Assert.hh, and avoids
 * inlining of exception throwing code.  Moreover, it also allows us to get \c
 * __FILE__ into play, which we couldn't do before because there was no way to
 * do string concatenation with it, and you couldn't do free store management
 * effectively in an inlined macro.  By going out to a global function, we pick
 * up the ability to formulate a more complete picture of the error, and provide
 * some optimization capability (both compile time and run time).
 *
 * Note also that at this juncture, we go ahead and drop all support for
 * compilers which are incapable of compiling exception code.  From here
 * forward, exceptions are assumed to be available.  And thereby we do our part
 * to promote "standard C++".
 *
 * \sa http://akrzemi1.wordpress.com/2013/01/04/preconditions-part-i/
 */
//===========================================================================//

class DLL_PUBLIC_dsxx assertion : public std::logic_error {
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
  explicit assertion(std::string const &msg)
      : std::logic_error(msg) { /* empty */
  }

  /*!
   * \brief Specialized constructor for rtt_dsxx::assertion class.
   *
   * This constructor creates a ds++ exception object.  This object is derived
   * form std::runtime_error and has identical functionality.  This constructor
   * is specialized for use by Draco DbC commands (Require, Ensure, Check, and
   * Insist).  It forms the error message from the test condition and the file
   * and line number of the DbC command.
   *
   * \param cond The expression that failed a DbC test.
   * \param file The source code file name that contains the DbC test.
   * \param line The source code line number that contains the DbC test.
   *
   * \sa \ref Draco_DBC, --with-dbc[=level], Require, Ensure, Check, Insist
   */
  assertion(std::string const &cond, std::string const &file, int const line)
      : std::logic_error(build_message(cond, file, line)) { /* empty */
  }

  /*! \brief Destructor for ds++/assertion class.
   * We do not allow the destructor to throw! */
  virtual ~assertion() throw();

  /*! Helper function to build error message that includes source file name and
   *  line number. */
  static std::string build_message(std::string const &cond,
                                   std::string const &file, int const line);
};

//---------------------------------------------------------------------------//
// FREE NAMESPACE FUNCTIONS
//---------------------------------------------------------------------------//

//! Throw a rtt_dsxx::assertion for Require, Check, Ensure.
DLL_PUBLIC_dsxx void toss_cookies(std::string const &cond,
                                  std::string const &file, int const line);

DLL_PUBLIC_dsxx void toss_cookies_ptr(char const *const cond,
                                      char const *const file, int const line);

//! Throw a rtt_dsxx::assertion if condition fails
DLL_PUBLIC_dsxx void check_cookies(bool cond, char const *cond_text,
                                   char const *file, int line);

//! Print error w/o throw
DLL_PUBLIC_dsxx void show_cookies(std::string const &cond,
                                  std::string const &file, int const line);
//! Throw a rtt_dsxx::assertion for Insist.
DLL_PUBLIC_dsxx void insist(std::string const &cond, std::string const &msg,
                            std::string const &file, int const line);

//! Pointer version of insist
DLL_PUBLIC_dsxx void insist_ptr(char const *const cond, char const *const msg,
                                char const *const file, int const line);

//! Check version of insist
DLL_PUBLIC_dsxx void check_insist(bool cond, char const *const condstr,
                                  std::string const &msg,
                                  char const *const file, int const line);

//! Check Pointer version of insist
DLL_PUBLIC_dsxx void check_insist_ptr(bool cond, char const *const condstr,
                                      char const *const msg,
                                      char const *const file, int const line);

//! Add hostname and pid to error messages.
DLL_PUBLIC_dsxx std::string verbose_error(std::string const &message);

} // end of rtt_dsxx

//---------------------------------------------------------------------------//
/*!
 * \page Draco_DBC Using the Draco Design-by-Contract Macros
 *
 * \section ddbc Using the Draco Design-by-Contract Macros
 *
 * The assertion macros are intended to be used for validating preconditions
 * which must be true in order for following code to be correct, etc.  For
 * example,
 *
 * \code
 * Assert( x > 0. );
 * y = sqrt(x);
 * \endcode
 *
 * If the assertion fails, the code should just bomb.  Philosophically, it
 * should be used to feret out bugs in preceding code, making sure that prior
 * results are within reasonable bounds before proceeding to use those results
 * in further computation, etc.
 *
 * These macros are provided to support the Design By Contract formalism.  The
 * activation of each macro is keyed off a bit in the DBC macro which can be
 * specified on the command line:
 *
 \verbatim
 Bit     DBC macro affected
 ---     ------------------
 0      Require
 1      Check
 2      Ensure, Remember
 3      (nothrow option)
 \endverbatim
 *
 * So for instance, \c -DDBC=7 turns them all on, \c -DDBC=0 turns them all off,
 * and \c -DDBC=1 turns on \c Require but turns off \c Check and \c Ensure.  The
 * default is to have them all enabled.  The fourth bit (3) can be toggled on to
 * disable the C++ thrown exception while keeping all of the DBC checks and
 * messages active.
 *
 * The \c Insist macro is akin to the \c Assert macro, but it provides the
 * opportunity to specify an instructive message.  The idea here is that you
 * should use Insist for checking things which are more or less under user
 * control.  If the user makes a poor choice, we "insist" that it be corrected,
 * providing a corrective hint.
 *
 * \note We provide a way to eliminate assertions, but not insistings.  The idea
 * is that \c Assert is used to perform sanity checks during program
 * development, which you might want to eliminate during production runs for
 * performance sake.  Insist is used for things which really really must be
 * true, such as "the file must've been opened", etc.  So, use \c Assert for
 * things which you want taken out of production codes (like, the check might
 * inhibit inlining or something like that), but use Insist for those things you
 * want checked even in a production code.
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
/*!
 * \def Insist_ptr(condition, message)
 *
 * Same as Insist, except that it uses char pointers, rather than strings.  This
 * is more efficient when inlined.
 */
//---------------------------------------------------------------------------//

// clang-format off

#if !defined(DBC)
#define DBC 7
#endif

//---------------------------------------------------------------------------//
// No-throw versions of DBC
//---------------------------------------------------------------------------//
#if DBC & 8

#if DBC & 1
#define REQUIRE_ON
#define Require(c) if (!(c)) rtt_dsxx::show_cookies( #c, __FILE__, __LINE__ )
#else
#define Require(c)
#endif

#if DBC & 2
#define CHECK_ON
#define Check(c) if (!(c)) rtt_dsxx::show_cookies( #c, __FILE__, __LINE__ )
#define Assert(c) if (!(c)) rtt_dsxx::show_cookies( #c, __FILE__, __LINE__ )
#else
#define Check(c)
#define Assert(c)
#endif

#if DBC & 4
#define ENSURE_ON
#define Ensure(c) if (!(c)) rtt_dsxx::show_cookies( #c, __FILE__, __LINE__ )
#else
#define Ensure(c)
#endif

//---------------------------------------------------------------------------//
// Always on
//---------------------------------------------------------------------------//
#define Insist(c,m) if (!(c)) rtt_dsxx::insist( #c, m, __FILE__, __LINE__ )
#define Insist_ptr(c,m) if (!(c)) rtt_dsxx::insist_ptr( #c, m, __FILE__, __LINE__ )

#elif DBC & 16

//---------------------------------------------------------------------------//
// Regular (exception throwing) versions of DBC, but with the check deferred.
// This eliminates numerous untestable branches for coverage analysis, but can
// be costly in run time.
// ---------------------------------------------------------------------------//

// Note that the !!(c) is needed to force the argument to a bool for certain
// objects, notably instances of the stream classes, which provide an operator!
// but no direct conversion to bool.

// Also note that Bad_Case is always off in this selection.

#if DBC & 1
#define REQUIRE_ON
#define Require(c) rtt_dsxx::check_cookies(!!(c), #c, __FILE__, __LINE__)
#else
#define Require(c)
#endif

#if DBC & 2
#define CHECK_ON
#define Check(c) rtt_dsxx::check_cookies(!!(c), #c, __FILE__, __LINE__)
#define Assert(c) rtt_dsxx::check_cookies(!!(c), #c, __FILE__, __LINE__)
#define Bad_Case(m)
#else
#define Check(c)
#define Assert(c)
#define Bad_Case(m)
#endif

#if DBC & 4
#define ENSURE_ON
#define Ensure(c) rtt_dsxx::check_cookies(!!(c), #c, __FILE__, __LINE__)
#else
#define Ensure(c)
#endif

//---------------------------------------------------------------------------//
// Always on
//---------------------------------------------------------------------------//
#define Insist(c, m) rtt_dsxx::check_insist(!!(c), #c, m, __FILE__, __LINE__)
#define Insist_ptr(c,m) rtt_dsxx::check_insist_ptr( !!(c), #c, m, __FILE__, __LINE__ )

#else // not DBC & 8 or DBC & 16

//---------------------------------------------------------------------------//
// Regular (exception throwing) versions of DBC
//---------------------------------------------------------------------------//
#if DBC & 1
#define REQUIRE_ON
#define Require(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ )
#else
#define Require(c)
#endif

#if DBC & 2
#define CHECK_ON
#define Check(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ )
#define Assert(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ )
#define Bad_Case(m) default: rtt_dsxx::toss_cookies( m, __FILE__, __LINE__ )
#else
#define Check(c)
#define Assert(c)
#define Bad_Case(m)
#endif

#if DBC & 4
#define ENSURE_ON
#define Ensure(c) if (!(c)) rtt_dsxx::toss_cookies( #c, __FILE__, __LINE__ )
#else
#define Ensure(c)
#endif

//---------------------------------------------------------------------------//
// Always on
//---------------------------------------------------------------------------//
#define Insist(c,m) if (!(c)) rtt_dsxx::insist( #c, m, __FILE__, __LINE__ )
#define Insist_ptr(c,m) if (!(c)) rtt_dsxx::insist_ptr( #c, m, __FILE__, __LINE__ )

#endif // DBC & 8

//---------------------------------------------------------------------------//
// If any of DBC is on, then make the remember macro active
//---------------------------------------------------------------------------//
#if DBC
#define REMEMBER_ON
#define Remember(c) c
#else
#define Remember(c)
#endif

#if defined(MSVC)
#pragma warning(pop)
#endif

// clang-format on

#endif // RTT_dsxx_Assert_HH

//---------------------------------------------------------------------------//
// end of ds++/Assert.hh
//---------------------------------------------------------------------------//
