//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/DBC_Ptr.hh
 * \author Paul Henning
 * \brief  Pointer-like class that is checked in DBC mode.
 * \note   Copyright &copy; 2005 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_DBC_Ptr_HH
#define rtt_dsxx_DBC_Ptr_HH

#include "Assert.hh"

/*!  
  \def DBC_Ptr
  This macro defines a pointer-like class for scalars.  The class is
  templated on the type, so you would make a declaration like:
  \code
  DBC_Ptr<Foo> foo_ptr(new Foo);
  \endcode 

  If DBC is on, DBC_Ptr will be the class rtt_dsxx::Safe_Ptr, which provides
  a variety of (possibly time and space expensive) diagnostics.  Otherwise,
  DBC_Ptr will be the class rtt_dsxx::Thin_Ptr, which should inline to being
  a bare pointer.
 
  \note 
    - do \em NOT point to arrays with these classes! The delete semantics
      are only for scalars.
    - these classes do \em NOT have automatic garabage collection!  You
      must call the \c delete_data() member functions to release memory.
    - like SP.hh, it is not hard to fool these classes... don't mix DBC_Ptr
      and rtt_dsxx::SP and/or raw pointer handles to the same data!
    - The definition and dependency graph in Doxygen-generated HTML is
      probably wrong!

  The intention of this class is to allow checked/safe development of 
  raw pointer like code.  When DBC is turned off for production code, all
  of the checking infrastructure disappears, and the compiler inlines the
  code back to raw pointers for efficiency.
*/

#if DBC

// ==========================================================================
//				  SAFE MODE
//			    Only when DBC is ON!
// ==========================================================================

#include "Safe_Ptr.hh"
#define DBC_Ptr rtt_dsxx::Safe_Ptr

#else

// ==========================================================================
//				  THIN MODE
//			    Only when DBC is OFF!
// ==========================================================================


#include "Thin_Ptr.hh"
#define DBC_Ptr rtt_dsxx::Thin_Ptr

#endif

#endif                          // rtt_dsxx_DBC_Ptr_HH

//---------------------------------------------------------------------------//
//                              end of ds++/DBC_Ptr.hh
//---------------------------------------------------------------------------//
