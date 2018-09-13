//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   quadrature/QIM_parser.cc
 * \author Kent Budge
 * \brief  Define a parse routine for quadrature interpolation model
 *         specifications.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "QIM.hh"

namespace rtt_quadrature {
using namespace rtt_parser;

//----------------------------------------------------------------------------//
/*!
 * /param tokens Token stream from which to parse the specification.
 *
 * /param qim Reference to a QIM into which the specification should be
 * stored. The routine checks that the QIM is set to END_QIM and reports a
 * semantic error if it is not. This simplifies checking for duplicate
 * specifications, by allowing the client to set the QIM to END_QIM before
 * beginning his parse.
 */
void parse_quadrature_interpolation_model(Token_Stream &tokens, QIM &qim) {
  tokens.check_semantics(qim == END_QIM,
                         "quadrature interpolation model already specified");

  Token token = tokens.shift();

  if (token.text() == "SN") {
    qim = SN;
  } else if (token.text() == "GQ1") {
    qim = GQ1;
  } else if (token.text() == "GQ2") {
    qim = GQ2;
  } else if (token.text() == "GQF") {
    qim = GQF;
  } else {
    tokens.check_semantics(false,
                           "unrecognized quadrature interpolation model");
  }
}

//----------------------------------------------------------------------------//
//! Provide a string representation of the provided quadrature enum.
std::string quadrature_interpolation_model_as_text(QIM q) {
  switch (q) {
  case SN:
    return "SN";
  case GQ1:
    return "GQ1";
  case GQ2:
    return "GQ2";
  case GQF:
    return "GQF";
  default:
    Insist(false, "bad case");
    return 0; // to kill warnings; never reached
  }
}

} // end namespace rtt_quadrature

//----------------------------------------------------------------------------//
// end of quadrature/QIM.hh
//----------------------------------------------------------------------------//
