//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/ExpInt.hh
 * \author Paul Talbot
 * \date   Tue Jul 26 14:48:13 MDT 2011
 * \brief  Declare the ExpInt function templates.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$ 
//---------------------------------------------------------------------------//
      

#ifndef special_functions_ExpInt_hh
#define special_functions_ExpInt_hh

namespace rtt_sf
{

//! Compute general exponential integral, order n, argument x \f$ E_n(x) \f$.
double En( unsigned const n, double const x);

//! Compute exponential integral, argument x \f$ Ei(x) \f$.
double Ei( double const x);

}

#endif //special_functions_ExpInt
