//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Gauss_Legendre.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright 2000-2010 Los Alamos National Security, LLC. All rights
 *         reserved. 
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include <numeric>

#include "Gauss_Legendre.hh"

#include "ds++/to_string.hh"
#include "gauleg.hh"

namespace rtt_quadrature
{
using namespace std;
using rtt_dsxx::to_string;

//---------------------------------------------------------------------------------------//
/* virtual */
string Gauss_Legendre::name() const
{
    return "Gauss-Legendre";
}

//---------------------------------------------------------------------------------------//
/* virtual */
string Gauss_Legendre::parse_name() const
{
    return "gauss legendre";
}

//---------------------------------------------------------------------------------------//
/* virtual */ string Gauss_Legendre::as_text(string const &indent) const
{
    string Result =
        indent + "type = gauss legendre" +
        indent + "  order = " + to_string(sn_order_) +
        indent + "end";

    return Result;
}

//---------------------------------------------------------------------------------------//
bool Gauss_Legendre::check_class_invariants() const
{
    return sn_order_>0 && sn_order_%2==0;
}

//---------------------------------------------------------------------------------------//
/* virtual */
vector<Ordinate>
Gauss_Legendre::create_level_ordinates_(double const norm) const
{
    // Preconditions checked in create_ordinate_set
    
    using rtt_dsxx::soft_equiv;

    double const mu1 = -1; // range of direction
    double const mu2 = 1;
    vector<double> mu, wt;
    unsigned const N = sn_order_;
    mu.reserve(N);
    wt.reserve(N);
    gauleg( mu1, mu2, mu, wt, sn_order_ );
    
    double sumwt = accumulate( wt.begin(),
                               wt.end(),
                               0.0 );  

    // If norm != 2.0 then renormalize the weights to the required values. 
    if( !soft_equiv(norm,2.0) ) 
    {
	double c = norm/sumwt;
	for ( size_t i=0; i < sn_order_; ++i )
	    wt[i] = c * wt[i];
    }
    
    // build the set of ordinates
    vector<Ordinate> Result( sn_order_ );
    for ( size_t i=0; i<sn_order_; ++i )
    {
	// This is a 1D set.
	Result[i] = Ordinate(mu[i], wt[i]);
    }

    return Result;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                       end of quadrature/Quadrature.cc
//---------------------------------------------------------------------------------------//
