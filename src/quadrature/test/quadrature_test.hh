//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/quadrature_test.hh
 * \author Kent G. Budge
 * \brief  Define class quadrature_test
 * \note   Copyright (C) 2007-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_quadrature_test_hh
#define quadrature_quadrature_test_hh

#include "ds++/UnitTest.hh"
#include "../Quadrature.hh"

namespace rtt_quadrature
{
using rtt_dsxx::UnitTest;

void quadrature_test(UnitTest &ut,
                     Quadrature &quadrature);

void quadrature_integration_test(UnitTest &ut,
                                 Quadrature &quadrature);


} // end namespace rtt_quadrature

#endif // quadrature_quadrature_test_hh

//---------------------------------------------------------------------------//
// end of quadrature/quadrature_test.hh
//---------------------------------------------------------------------------//
