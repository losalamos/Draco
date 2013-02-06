//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Double_Gauss.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval dobule Gauss-Legendre quadrature set.
 * \note   Copyright 2000-2010 Los Alamos National Security, LLC. All rights
 *         reserved. 
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include <numeric>

#include "Double_Gauss.hh"

#include "ds++/to_string.hh"
#include "gauleg.hh"

namespace rtt_quadrature
{
using namespace std;
using rtt_dsxx::to_string;

//---------------------------------------------------------------------------------------//
/* virtual */
string Double_Gauss::name() const
{
    return "Double-Gauss";
}

//---------------------------------------------------------------------------------------//
/* virtual */
string Double_Gauss::parse_name() const
{
    return "double gauss";
}

//---------------------------------------------------------------------------------------//
/* virtual */
unsigned Double_Gauss::number_of_levels() const
{
    return sn_order();
}

//---------------------------------------------------------------------------------------//
/* virtual */ string Double_Gauss::as_text(string const &indent) const
{
    string Result =
        indent + "type = double gauss" +
        indent + "  order = " + to_string(sn_order()) +
        indent + "end";

    return Result;
}

//---------------------------------------------------------------------------------------//
bool Double_Gauss::check_class_invariants() const
{
    return sn_order()>0 && sn_order()%2==0;
}

//---------------------------------------------------------------------------------------//
/* virtual */
vector<Ordinate>
Double_Gauss::create_level_ordinates_(double const norm) const
{
    // Preconditions checked in create_ordinate_set

    unsigned const numGaussPoints = sn_order();
    unsigned const n(numGaussPoints);
    unsigned const n2(n/2);

    // size the data vectors

    vector<double> mu(numGaussPoints);
    vector<double> wt(numGaussPoints);

    if (n2 == 1)
    {
        // 2-point double Gauss is just Gauss

        double const mu1(-1); // range of direction
        double const mu2(1);
        gauleg( mu1, mu2, mu, wt, numGaussPoints );

    }
    else
    {
        // Create an N/2-point Gauss quadrature on [-1,1]

        Check( n2%2 == 0 );

        double const mu1(-1); // range of direction
        double const mu2(1);
        std::vector< double > muH;
        std::vector< double > wtH;
        gauleg( mu1, mu2, muH, wtH, n2 );

        // map the quadrature onto the two half-ranges

        for (unsigned m=0; m<n2; ++m)
        {
            // Map onto [-1,0] then skew-symmetrize
            // (ensuring ascending order on [-1, 1])
            
    
            mu[m] = 0.5*(muH[m] - 1.0);
            wt[m] = 0.5*wtH[m];

            mu[n-m-1] = -mu[m];
            wt[n-m-1] =  wt[m];
        }
    }

    // Compute and store the sum of the weights

    double sumwt( std::accumulate( wt.begin(), wt.end(), // range
                                   0.0 ) );              // init value.

    // Sanity Checks: always none at present

    if( !soft_equiv(norm,2.0) )
    {
        double c = norm/sumwt;
        for ( size_t i=0; i < numGaussPoints; ++i )
            wt[i] = c * wt[i];
    }

    // build the set of ordinates
    vector<Ordinate> Result( numGaussPoints ); // sn_order
    for ( size_t i=0; i<numGaussPoints; ++i )
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
