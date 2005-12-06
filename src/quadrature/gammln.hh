//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/gammln.hh
 * \author Kent Budge
 * \date   Tue Sep 14 13:52:58 2004
 * \brief  Log of gamma function
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef utils_gammln_hh
#define utils_gammln_hh

namespace rtt_utils
{

//! Log of gamma function
template <class Field> Field gammln(const Field xx)
{
  Field x,y,tmp,ser;
  int j;

  const double cof[6] = {76.18009172947146, -86.50532032941677,
			 24.01409824083091, -1.2317395572450155,
			 0.12086500973866179e-2, -0.5395239384953e-5};
  
  y = x = xx;
  tmp = x+5.5;
  tmp -= (x+0.5)*std::log(tmp);
  ser = 1.000000000190015;
  for (j=0; j<=5; j++) y+=1.0, ser += cof[j]/y;
  return -tmp + std::log(2.5066282746310005*ser/x);
}

} // end namespace rtt_utils

#endif // utils_gammln_hh

//---------------------------------------------------------------------------//
//              end of utils/gammln.hh
//---------------------------------------------------------------------------//
