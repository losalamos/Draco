//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/to_string.hh
 * \author Kent Budge
 * \brief  Define class to_string
 * \note   Copyright (C) 2007 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_to_string_hh
#define dsxx_to_string_hh

#include <sstream>
#include <string>

namespace rtt_dsxx
{
using std::string;

//---------------------------------------------------------------------------//
// Simple function which converts a number into a string.
// http://public.research.att.com/~bs/bs_faq2.html
//---------------------------------------------------------------------------//

template< class T >
string to_string(T num, unsigned precision = 0 )
{
    using namespace std;
    
    std::stringstream s;
    if (precision) s.precision(precision);
    s << num ;
    return s.str();
}

} // end namespace rtt_dsxx

#endif // dsxx_to_string_hh

//---------------------------------------------------------------------------//
//              end of ds++/to_string.hh
//---------------------------------------------------------------------------//
