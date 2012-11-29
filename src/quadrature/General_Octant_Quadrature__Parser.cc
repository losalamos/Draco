//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/General_Octant_Quadrature.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright 2000-2010 Los Alamos National Security, LLC. All rights
 *         reserved. 
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include "parser/utilities.hh"
#include "General_Octant_Quadrature.hh"

namespace rtt_quadrature
{
using namespace rtt_parser;

//---------------------------------------------------------------------------------------//
/*static*/
SP<Quadrature> General_Octant_Quadrature::parse(Token_Stream &tokens)
{
    Token token = tokens.shift();
    tokens.check_syntax(token.text()=="number of ordinates",
                        "expected number of ordinates");

    unsigned N = parse_positive_integer(tokens);

    token = tokens.shift();
    tokens.check_syntax(token.text()=="number of levels",
                        "expected number of levels");

    unsigned number_of_levels = parse_unsigned_integer(tokens);

    token = tokens.lookahead();
    
    QIM qim = SN;
    if (token.text()=="interpolation algorithm")
    {
        tokens.shift();
        token = tokens.shift();
        if (token.text()=="SN")
        {
            // default
        }
        else if (token.text()=="GALERKIN")
        {
            qim = GQ;
        }
    }

    vector<double> mu(N), eta(N), xi(N), wt(N);

    for (unsigned i=0; i<N; ++i)
    {
        mu[i] = parse_real(tokens);
        eta[i] = parse_real(tokens);
        xi[i] = parse_real(tokens);
        wt[i] = parse_real(tokens);
    }

    tokens.check_syntax(tokens.shift().type()==END, "missing end?");

    return SP<Quadrature>(new General_Octant_Quadrature(mu, eta, xi, wt,
                                                        number_of_levels,
                                                        qim));
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                       end of quadrature/Quadrature.cc
//---------------------------------------------------------------------------------------//
