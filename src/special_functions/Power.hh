//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   special_functions/Power.hh
 * \author Mike Buksas
 * \date   Thu Jul 20 17:23:31 2006
 * \brief  
 * \note   Copyright © 2006 Los Alamos National Security, LLC

 Use meta-programming to generate an efficient routine to compute integer
 powers. 

 E.g.: Power<4>(x) computes the fourth power of x. The function constructed at
 compile time should be equivalent to the following:

 v1 = x * x
 v2 = v1 * v1;
 return v2;

 Assuming that the compiler is inlining aggressively.

 Likewise, Power<7>(x) should generate code equivalent to:

 v1 = x * x
 v2 = v1 * x;
 v3 = v1 * v2;
 return v3;

 The meta-algorithm is based on the Russian Peasant Algorithm, modified to be
 recursive, since this is required for template meta-programming in C++.
 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef special_functions_Power_hh
#define special_functions_Power_hh

namespace rtt_sf
{

namespace {

/* Protect the implementation detail of struct P from accidental usage
 * outside of this file.
 */

template <int N, typename F>
struct P
{
    static F compute(F x, F p)
    {
        x *= x;
        if ((N/2)*2 == N) return P<N/2, F>::compute(x, p);
        else              return P<N/2, F>::compute(x, x*p);
    }
};

template <typename F>
struct P<0,F>
{
    static F compute(F x, F p) { return p; }
};

}
          

template <int N, typename F>
F Power(F x) {
    if ((N/2)*2 == N) return Power<N/2>(x*x);
    else              return P<N/2,F>::compute(x,x);
}


} // end namespace rtt_sf

#endif // special_functions_Power_hh

//---------------------------------------------------------------------------//
//              end of special_functions/Power.hh
//---------------------------------------------------------------------------//
