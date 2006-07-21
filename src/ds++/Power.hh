//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Power.hh
 * \author Mike Buksas
 * \date   Thu Jul 20 17:23:31 2006
 * \brief  
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_Power_hh
#define dsxx_Power_hh

namespace rtt_dsxx
{

template <int N, typename F>
struct RP2
{

    static double compute(F x, F p)
    {

        x *= x;
        if ((N/2)*2 == N) return RP2<N/2, F>::compute(x, p);
        else              return RP2<N/2, F>::compute(x, x*p);

    }

};


template <int N, typename F>
struct RP
{

    static double compute(F x)
    {

        if ((N/2)*2 == N) return RP <N/2, F>::compute(x*x);
        else              return RP2<N/2, F>::compute(x,x);

    }

};

template <>
template <typename F>
struct RP2<1, F>
{
    static double compute(F x, F p) { return x*x*p; }
};
    
template <>
template <typename F>
struct RP2<0,F>
{
    static double compute(F x, F p) { return p; }
};
          


template <int N, typename F>
F Power(F x) { return RP<N, F>::compute(x); }


} // end namespace rtt_dsxx

#endif // dsxx_Power_hh

//---------------------------------------------------------------------------//
//              end of ds++/Power.hh
//---------------------------------------------------------------------------//
