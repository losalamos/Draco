//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   shared_lib/test/Foo.hh
 * \author Rob Lowrie
 * \date   Wed Dec 29 11:24:17 2004
 * \brief  Header for Foo.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#ifndef shared_lib_test_Foo_hh
#define shared_lib_test_Foo_hh

#include "Foo_Base.hh"

namespace rtt_shared_lib_test {

//============================================================================//
/*!
  \class Foo
  \brief A simple class that may be defined in a shared library.

  When loaded as a shared lib, objects of this class are created and
  destroyed via the functions defined in creator_destroyer.cc.
  
  For more information, see the "C++ dlopen mini HOWTO" by Aaron Isotton.
 */
//============================================================================//

class Foo : public Foo_Base {
private:
  // DATA

  // A base value.
  double d_base;

public:
  // CREATORS

  // Constructor.
  Foo(const double x);

  // Compute something.
  double compute(const double x) const;
};

} // end namespace rtt_shared_lib_test

#endif // shared_lib_test_foo_hh

//----------------------------------------------------------------------------//
//              end of shared_lib/test/Foo.hh
//----------------------------------------------------------------------------//
