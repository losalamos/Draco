//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   shared_lib/test/Foo_Base.hh
 * \author Rob Lowrie
 * \date   Wed Dec 29 11:24:17 2004
 * \brief  Header for Foo_Base.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#ifndef shared_lib_test_Foo_Base_hh
#define shared_lib_test_Foo_Base_hh

namespace rtt_shared_lib_test {

//============================================================================//
/*!
  \class Foo_Base
  \brief An abstract base class.

  This class may act as an abstract base class for derived classes that are
  defined within a shared library.
  
  For more information, see the "C++ dlopen mini HOWTO" by Aaron Isotton.
 */
//============================================================================//

class Foo_Base {
public:
  // Must define a virual destructor.  It's implementation should be linked
  // in with the executable.
  virtual ~Foo_Base() = 0;

  // Compute something.
  virtual double compute(const double x) const = 0;
};

} // end namespace rtt_shared_lib_test

#endif // shared_lib_test_foo_hh

//----------------------------------------------------------------------------//
//              end of shared_lib/test/Foo_Base.hh
//----------------------------------------------------------------------------//
