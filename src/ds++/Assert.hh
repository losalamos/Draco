//----------------------------------*-C++-*----------------------------------//
// Assert.hh
// Geoffrey Furnish
// 6 December 1993
//---------------------------------------------------------------------------//
// @> Header file for Assert macro, and friends.
//---------------------------------------------------------------------------//
// This is my own personal shot at doing Assert right.  The idea for this was
// gleaned from a discussion with Mark Gray after he read something in
// "Writing Solid Code" about the assert macro being busted.  Should really
// go look at that discussion in depth at some point to see if I really need
// this, and if I should embellish it somehow.
//---------------------------------------------------------------------------//

#ifndef __ds_Assert_hh__
#define __ds_Assert_hh__

//---------------------------------------------------------------------------//
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
//---------------------------------------------------------------------------//

#include "config.hh"

NAMESPACE_DS_BEG

//===========================================================================//
// class assertion - exception notification class for assertions

// This class should really be derived from std::runtime_error, but
// unfortunately we don't have good implementation of the library standard
// yet, on compilers other than KCC.  So, this class will keep with the
// "what" method evidenced in the standard, but dispense with inheriting from
// classes for which we don't have implementations...
//===========================================================================//

class assertion /* : std::runtime_error */
{
    char *msg;
  public:
    assertion( const char *cond, const char *file, int line );
    assertion( const char *m );
    ~assertion() { delete[] msg; }
    const char* what() const { return msg; }
};

void toss_cookies( const char *cond, const char *file, int line );
void insist( const char *cond, const char *msg, const char *file, int line );

NAMESPACE_DS_END

//---------------------------------------------------------------------------//
// The Assert macro is intended to be used for validating preconditions which
// must be true in order for following code to be correct, etc.  For example,
// Assert( x > 0. ); y = sqrt(x);  If the assertion fails, the code should
// just bomb.  Philosophically, it should be used to feret out bugs in
// preceding code, making sure that prior results are within reasonable
// bounds before proceeding to use those results in further computation, etc.
//---------------------------------------------------------------------------//

#ifdef NOASSERT
#define Assert(c) 
#else
#define Assert(c) if (!(c)) dsxx::toss_cookies( #c, __FILE__, __LINE__ );
#endif

//---------------------------------------------------------------------------//
// The Insist macro is akin to the Assert macro, but it provides the
// opportunity to specify an instructive message.  The idea here is that you
// should use Insist for checking things which are more or less under user
// control.  If the user makes a poor choice, we "insist" that it be
// corrected, providing a corrective hint.
//---------------------------------------------------------------------------//

#define Insist(c,m) if (!(c)) dsxx::insist( #c, m, __FILE__, __LINE__ );

//---------------------------------------------------------------------------//
// NOTE:  We provide a way to eliminate assertions, but not insistings.  The
// idea is that Assert is used to perform sanity checks during program
// development, which you might want to eliminate during production runs for
// performance sake.  Insist is used for things which really really must be
// true, such as "the file must've been opened", etc.  So, use Assert for
// things which you want taken out of production codes (like, the check might
// inhibit inlining or something like that), but use Insist for those things
// you want checked even in a production code.
//---------------------------------------------------------------------------//

#endif				// __ds_Assert_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/Assert.hh
//---------------------------------------------------------------------------//
