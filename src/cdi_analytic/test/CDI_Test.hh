//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/CDI_Test.hh
 * \author Thomas M. Evans
 * \date   Wed Aug 29 11:16:17 2001
 * \brief  Test services for cdi_analytic package.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_test_CDI_Test_hh__
#define __cdi_analytic_test_CDI_Test_hh__

#include "../Analytic_Models.hh"

namespace rtt_cdi_analytic_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

bool fail(int line);

bool fail(int line, char *file);

#define FAILURE fail(__LINE__, __FILE__);

//===========================================================================//
// USER-DEFINED ANALYTIC_OPACITY_MODEL
//===========================================================================//

class Marshak_Model : public rtt_cdi_analytic::Analytic_Opacity_Model
{
    double a;
  public:
    Marshak_Model(double a_) : a(a_) {/*...*/}

    double calculate_opacity(double T, double rho) const
    {
	return a / (T * T * T);
    }
};

} // end namespace rtt_cdi_analytic_test

#endif                          // __cdi_analytic_test_CDI_Test_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_analytic/test/CDI_Test.hh
//---------------------------------------------------------------------------//
