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

#include <iostream.h>
#include <stdlib.h>

//---------------------------------------------------------------------------//
// The Assert macro is intended to be used for validating preconditions which
// must be true in order for following code to be correct, etc.  For example,
// Assert( x > 0. ); y = sqrt(x);  If the assertion fails, the code should
// just bomb.  Philosophically, it should be used to feret out bugs in
// preceding code, making sure that prior results are within reasonable
// bounds before proceeding to use those results in further computation, etc.
//---------------------------------------------------------------------------//

#ifdef NOASSERT
#define Assert(a) 
#else
#define Assert(a) if (!(a)) \
{ cerr << "Assertion " << #a << " failed in " << __FILE__ \
       << " at line " << __LINE__ << endl << flush; \
throw( "Assertion " #a " failed in " __FILE__ ); }
#endif

//---------------------------------------------------------------------------//
// The Insist macro is akin to the Assert macro, but it provides the
// opportunity to specify an instructive message.  The idea here is that you
// should use Insist for checking things which are more or less under user
// control.  If the user makes a poor choice, we "insist" that it be
// corrected, providing a corrective hint.
//---------------------------------------------------------------------------//

#define Insist(a,b) if (!(a)) \
{ cerr << "Insisting that: " << #a << " , in " << __FILE__ \
       << " at line " << __LINE__ << endl << flush; \
throw( b ); }

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

//---------------------------------------------------------------------------//
// The DIE macro is for putting into function bodies which are only place
// holders until they can be actually implemented.  This will encourage that
// process :-).
//---------------------------------------------------------------------------//

#define DIE \
{ cerr << "Function containing line " << __LINE__ << " in file " \
       << __FILE__ << " not implemented yet.\n"; \
  abort(); }

// The Tcl_Assert macro is pretty much like the Insist macro, except that it
// is intended to be used with Tcl.

#define Tcl_Assert(a,b) if (!(a)) { interp->result = b; return TCL_ERROR; }

//---------------------------------------------------------------------------//
// Okay, now we try to do something with exceptions.  The idea here is that
// an exception is something which had better not ever get compiled out.  The
// above asserions could conceivably get compiled out to nothing if the user
// wants to play games with CPP.  But an exception is a language
// feature--it's supposed to work no matter what.

// Alas, at the time of this writing, June 1994, most commercially available
// compilers on machines I am using, don't have exceptions.  So for the time
// being, I have to play CPP games here too.  Worse, both "throw" and "catch"
// have semantics that are beyond CPP's poor pea brain.  So I have to split
// it into throw/Throw and catch/Catch to handle the different syntactic
// forms.  Throw is for "re-throwing", and Catch is for catching with a
// variable spec, which must be defined in order for the following code to
// compile.  Note however, that the Catch variable is declared at one scope
// level higher than the catching code, which is different from how real C++
// works, so that variable name collisions are possible.  grrrr.  The only
// good solution to all this is to use exception capable compilers
// exclusively.  Some day ...
//---------------------------------------------------------------------------//

// Okay, HP C++ has exceptions, but most other compilers don't.  Will
// embellish this as necessary.

#ifdef __hpux
#if defined(__GNUC__) || defined(__lucid) || defined(__CENTERLINE__) \
|| defined(CENTERLINE_CLPP)
#define NO_XCPT
#endif
#else
// #ifndef __linux
#define NO_XCPT
// #endif
#endif

// Kuck and Associates C++ 3.0+ has exceptions on all platforms.

#ifdef __KCC
#undef NO_XCPT
#endif

#ifdef _POWER
#undef NO_XCPT
#endif

#ifdef __DECCXX
#undef NO_XCPT
#endif

#ifdef HAVE_XCPT
#undef NO_XCPT
#endif

#ifdef NO_XCPT
#define try
#define throw(a) \
{ cerr << "THROW: " << #a << " from " << __FILE__ \
       << " line " << __LINE__ << endl << flush; abort(); }
#define catch(a) if (0)
#define Catch(a) a; if (0)
#define Throw
#else
#define Catch(a) catch(a)
#define Throw throw
#endif

#endif				// __ds_Assert_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/Assert.hh
//---------------------------------------------------------------------------//
