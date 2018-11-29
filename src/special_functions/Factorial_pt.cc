//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/factorial_pt.cc
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide explicit instantiations of templatized factorial function. 
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Factorial.i.hh"

namespace rtt_sf {

//---------------------------------------------------------------------------//
// Make factorial valid only for int and unsigned.

template DLL_PUBLIC_special_functions unsigned factorial(unsigned const k);

template DLL_PUBLIC_special_functions int factorial(int const k);

template DLL_PUBLIC_special_functions long factorial(long const k);

template DLL_PUBLIC_special_functions double
factorial_fraction(unsigned const k, unsigned const l);

template DLL_PUBLIC_special_functions double factorial_fraction(int const k,
                                                                int const l);

template DLL_PUBLIC_special_functions double factorial_fraction(long const k,
                                                                long const l);

} // end namespace rtt_sf

//---------------------------------------------------------------------------//
// end of sf/factorial_pt.cc
//---------------------------------------------------------------------------//
