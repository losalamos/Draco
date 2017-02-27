//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/to_string.hh
 * \author Kent Budge
 * \brief  Define class to_string
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_to_string_hh
#define dsxx_to_string_hh

#include <sstream>
#include <string>

namespace rtt_dsxx {
using std::string;

//---------------------------------------------------------------------------//
// Simple function which converts a number into a string.
// http://public.research.att.com/~bs/bs_faq2.html
//---------------------------------------------------------------------------//

template <class T>
string to_string(T const num, unsigned int const precision = 23) {
  std::stringstream s;
  s.precision(precision);
  s << num;
  return s.str();
}

} // end namespace rtt_dsxx

#endif // dsxx_to_string_hh

//---------------------------------------------------------------------------//
// end of ds++/to_string.hh
//---------------------------------------------------------------------------//
