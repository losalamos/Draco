//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/cdi_analytic_test.hh
 * \author Thomas M. Evans
 * \date   Mon Sep 24 12:04:00 2001
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_test_hh__
#define __cdi_analytic_test_hh__

#include "../Analytic_Models.hh"
#include "ds++/Packing_Utils.hh"
#include <iostream>
#include <vector>

namespace rtt_cdi_analytic_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_cdi_analytic_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_cdi_analytic_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain draco tests except that the output must be of
// the form "Test: pass/fail"

bool fail(int line);

bool fail(int line, char *file);

bool pass_msg(const std::string &);

bool fail_msg(const std::string &);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

extern bool passed;

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

    double calculate_opacity(double T, double rho, double nu) const
    {
        return calculate_opacity(T, rho);
    }

    std::vector<double> get_parameters() const
    {
	return std::vector<double>(1, a);
    }

    std::vector<char> pack() const
    {
	rtt_dsxx::Packer packer;
	std::vector<char> p(sizeof(double) + sizeof(int));
	packer.set_buffer(p.size(), &p[0]);
	int indicator = 10;
	packer << indicator << a;
	return p;
    }
};

} // end namespace rtt_cdi_analytic_test

#define ITFAILS    rtt_cdi_analytic_test::fail(__LINE__);
#define FAILURE    rtt_cdi_analytic_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_cdi_analytic_test::pass_msg(a);
#define FAILMSG(a) rtt_cdi_analytic_test::fail_msg(a);

#endif                          // __cdi_analytic_test_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_analytic/test/cdi_analytic_test.hh
//---------------------------------------------------------------------------//
