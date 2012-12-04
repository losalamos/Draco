//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/QIM_parser.cc
 * \author Kent Budge
 * \brief  Define a parse routine for quadrature interpolation model specifications.
 * \note   © Copyright 2006-2012 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------------------//
// $Id: QIM.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include "QIM.hh"

namespace rtt_quadrature
{
using namespace rtt_parser;

//---------------------------------------------------------------------------------------//
/*!
 * /param tokens Token stream from which to parse the specification.
 *
 * /param qim Reference to a QIM into which the specification should be
 * stored. The routine checks that the QIM is set to END_QIM and reports a
 * semantic error if it is not. This simplifies checking for duplicate
 * specifications, by allowing the client to set the QIM to END_QIM before
 * beginning his parse.
 */
void parse_quadrature_interpolation_model(Token_Stream &tokens,
                                          QIM &qim)
{
    tokens.check_semantics(qim==END_QIM,
                           "quadrature interpolation model already specified");
    
    Token token = tokens.shift();

    tokens.check_semantics(token.text()=="SN" || token.text()=="GALERKIN",
                           "unrecognized quadrature interpolation model");
    
    if (token.text()=="SN")
    {
        qim = SN;
    }
    else if (token.text()=="GALERKIN")
    {
        qim = GQ;
    }
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//              end of quadrature/QIM.hh
//---------------------------------------------------------------------------------------//
