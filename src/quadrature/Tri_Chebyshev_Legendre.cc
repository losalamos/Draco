//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Tri_Chebyshev_Legendre.cc
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------------------//
// $Id: Tri_Chebyshev_Legendre.cc 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

// #include <iostream>
// #include <iomanip>
// #include <cmath>
// #include <algorithm>

#include "Tri_Chebyshev_Legendre.hh"

#include "gauleg.hh"
#include "ds++/to_string.hh"

namespace rtt_quadrature
{
using namespace rtt_dsxx;

//---------------------------------------------------------------------------------------//
string Tri_Chebyshev_Legendre::name() const { return "Tri Chebyshev Legendre"; }

//---------------------------------------------------------------------------------------//
string Tri_Chebyshev_Legendre::parse_name()  const { return "tri cl"; }
    
//---------------------------------------------------------------------------------------//
Quadrature::Quadrature_Class Tri_Chebyshev_Legendre::quadrature_class() const
{
    return TRIANGLE_QUADRATURE;
}

//---------------------------------------------------------------------------------------//
unsigned Tri_Chebyshev_Legendre::number_of_levels() const { return sn_order_; }
    
//---------------------------------------------------------------------------------------//
string Tri_Chebyshev_Legendre::as_text(string const &indent) const
{
    string Result =
        indent + "type = tri cl" +
        indent + "  order = " + to_string(sn_order_) +
        Octant_Quadrature::as_text(indent);

    return Result;
}

//---------------------------------------------------------------------------------------//
void
Tri_Chebyshev_Legendre::create_octant_ordinates_(vector<double> &mu,
                                                 vector<double> &eta,
                                                 vector<double> &wt) const
{
    using rtt_dsxx::soft_equiv;

    // The number of quadrature levels is equal to the requested SN order.
    size_t levels = sn_order_;

    // We build the 3-D first, then edit as appropriate.

    size_t numOrdinates = levels*(levels+2)/8;

    // Force the direction vectors to be the correct length.
    mu.resize(numOrdinates);
    eta.resize(numOrdinates);
    wt.resize(numOrdinates);

    double const mu1 = -1; // range of direction
    double const mu2 = 1;
    vector<double> GLmu(levels);
    vector<double> GLwt(levels);
    gauleg( mu1, mu2, GLmu, GLwt, sn_order_ );

    // NOTE: this aligns the gauss points with the x-axis (r-axis in cylindrical coords)

    unsigned icount=0;

    for (unsigned i=0; i<levels/2; ++i)
    {
        double xmu=GLmu[i];
        double xwt=GLwt[i];
        double xsr=sqrt(1.0-xmu*xmu);

        unsigned const k=2*(i+1);
            
        for (unsigned j=0; j<k/2; ++j)
        {
            unsigned ordinate=icount;
                
            mu[ordinate] =  xsr*cos(rtt_units::PI*(2.0*j+1.0)/k/2.0);
            eta[ordinate] = xsr*sin(rtt_units::PI*(2.0*j+1.0)/k/2.0);
            wt[ordinate] = xwt/k;

            ++icount;
        }
    }
}


} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                 end of Tri_Chebyshev_Legendre.cc
//---------------------------------------------------------------------------------------//
