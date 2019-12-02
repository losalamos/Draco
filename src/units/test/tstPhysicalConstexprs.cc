//-----------------------------------*-C++-*----------------------------------//
/*!
 * \file   units/test/tstPhysicalConstexprs.cc
 * \author Ben R. Ryan <brryan@lanl.gov>
 * \date   2017 Feb 10
 * \brief  Implementation file for PhysicalConstexprs class test
 * \note   Copyright (C) 2017-2019 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstexprs.hh"
#include <iomanip>
#include <sstream>

namespace rtt_units_test {

//! Helper function for comparing PhysicalConstexprs values to expectation
void compare_constant(const std::string name, double expected, double value_a,
                      double value_b, rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using std::endl;
  using std::setprecision;

  if (soft_equiv(expected, value_a) && soft_equiv(expected, value_b)) {
    PASSMSG("Found expected value for " + name);
  } else {
    std::ostringstream msg;
    msg << "Did not find expected value for " + name + "." << endl
        << "\tThe value expected was " << setprecision(16) << expected << "."
        << endl
        << "\tThe values returned were " << setprecision(16) << value_a
        << " and " << setprecision(16) << value_b << "." << endl;
    FAILMSG(msg.str());
  }
}

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

//! Tests the SI unit system
void test_si(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;
  using std::endl;
  using std::ostringstream;

  constexpr rtt_units::PhysicalConstexprs<rtt_units::SI> pc;

  compare_constant("Avogadro's number", 6.02214076e23, pc.avogadro(), pc.Na(),
                   ut);
  compare_constant("fine structure constant", 0.0072973525693,
                   pc.fineStructure(), pc.alpha(), ut);
  compare_constant("pi", 3.141592653589793238462643383279, pc.pi(), pc.pi(),
                   ut);
  compare_constant("Planck constant", 6.62607015e-34, pc.planck(), pc.h(), ut);
  compare_constant("molar gas constant", 8.3144621618, pc.gasConstant(), pc.R(),
                   ut);
  compare_constant("Boltzmann constant", 1.602176797090563e-25, pc.boltzmann(),
                   pc.k(), ut);
  compare_constant("electron charge", 1.602176634e-19, pc.electronCharge(),
                   pc.e(), ut);
  compare_constant("speed of light", 299792458, pc.speedOfLight(), pc.c(), ut);
  compare_constant("Stefan-Boltzmann constant", 5.670374419e-08,
                   pc.stefanBoltzmann(), pc.sigma(), ut);
  compare_constant("Gravitational constant", 6.6743e-11,
                   pc.gravitationalConstant(), pc.G(), ut);
  compare_constant("Standard acceleration of gravity", 9.80665,
                   pc.accelerationFromGravity(), pc.g(), ut);
  compare_constant("Faraday constant", 96485.33212, pc.faradayConstant(),
                   pc.F(), ut);
  compare_constant("Permeability of vacuum", 1.25663706212e-06,
                   pc.permeabilityOfVacuum(), pc.mu0(), ut);
  compare_constant("Permittivity of vacuum", 8.854187812800001e-12,
                   pc.permittivityOfVacuum(), pc.epsi0(), ut);
  compare_constant("Classical electron radius", 2.8179403262e-13,
                   pc.classicalElectronRadius(), pc.re(), ut);
  compare_constant("Electron mass", 9.1093837015e-28, pc.electronMass(),
                   pc.me(), ut);
  compare_constant("Proton mass", 1.67262192369e-24, pc.protonMass(), pc.mp(),
                   ut);
  compare_constant("Electron volt", 1.602176634e-28, pc.electronVolt(), pc.eV(),
                   ut);
}

//! Tests the cgs unit system
void test_cgs(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;
  using std::endl;
  using std::ostringstream;

  constexpr rtt_units::PhysicalConstexprs<rtt_units::CGS> pc;

  compare_constant("Avogadro's number", 6.02214076e23, pc.avogadro(), pc.Na(),
                   ut);
  compare_constant("fine structure constant", 0.0072973525693,
                   pc.fineStructure(), pc.alpha(), ut);
  compare_constant("pi", 3.141592653589793238462643383279, pc.pi(), pc.pi(),
                   ut);
  compare_constant("Planck constant", 6.626070149999998e-27, pc.planck(),
                   pc.h(), ut);
  compare_constant("molar gas constant", 83144621.618, pc.gasConstant(), pc.R(),
                   ut);
  compare_constant("Boltzmann constant", 1.602176797090563e-25, pc.boltzmann(),
                   pc.k(), ut);
  compare_constant("electron charge", 4.803204712570263e-10,
                   pc.electronCharge(), pc.e(), ut);
  compare_constant("speed of light", 29979245800., pc.speedOfLight(), pc.c(),
                   ut);
  compare_constant("Stefan-Boltzmann constant", 5.670374419e-05,
                   pc.stefanBoltzmann(), pc.sigma(), ut);
  compare_constant("Gravitational constant", 6.6743e-08,
                   pc.gravitationalConstant(), pc.G(), ut);
  compare_constant("Standard acceleration of gravity", 980.665,
                   pc.accelerationFromGravity(), pc.g(), ut);
  compare_constant("Faraday constant", 8.667420510246794e+16,
                   pc.faradayConstant(), pc.F(), ut);
  compare_constant("Permeability of vacuum", 12.5663706212,
                   pc.permeabilityOfVacuum(), pc.mu0(), ut);
  compare_constant("Permittivity of vacuum", 0.07953848254861551,
                   pc.permittivityOfVacuum(), pc.epsi0(), ut);
  compare_constant("Classical electron radius", 2.8179403262e-13,
                   pc.classicalElectronRadius(), pc.re(), ut);
  compare_constant("Electron mass", 9.1093837015e-28, pc.electronMass(),
                   pc.me(), ut);
  compare_constant("Proton mass", 1.67262192369e-24, pc.protonMass(), pc.mp(),
                   ut);
  compare_constant("Electron volt", 1.602176634e-12, pc.electronVolt(), pc.eV(),
                   ut);
}

//! Tests the cgsh unit system
void test_cgsh(rtt_dsxx::UnitTest &ut) {
  using rtt_dsxx::soft_equiv;
  using rtt_units::PI;
  using std::endl;
  using std::ostringstream;

  constexpr rtt_units::PhysicalConstexprs<rtt_units::CGSH> pc;

  compare_constant("Avogadro's number", 6.02214076e23, pc.avogadro(), pc.Na(),
                   ut);
  compare_constant("fine structure constant", 0.0072973525693,
                   pc.fineStructure(), pc.alpha(), ut);
  compare_constant("pi", 3.141592653589793238462643383279, pc.pi(), pc.pi(),
                   ut);
  compare_constant("Planck constant", 6.626070149999998e-35, pc.planck(),
                   pc.h(), ut);
  compare_constant("molar gas constant", 0.09648533664909333, pc.gasConstant(),
                   pc.R(), ut);
  compare_constant("Boltzmann constant", 1.602176797090563e-25, pc.boltzmann(),
                   pc.k(), ut);
  compare_constant("electron charge", 4.803204712570263e-10,
                   pc.electronCharge(), pc.e(), ut);
  compare_constant("speed of light", 299.792458, pc.speedOfLight(), pc.c(), ut);
  compare_constant("Stefan-Boltzmann constant", 1.028301235680117,
                   pc.stefanBoltzmann(), pc.sigma(), ut);
  compare_constant("Gravitational constant", 6.6743e-24,
                   pc.gravitationalConstant(), pc.G(), ut);
  compare_constant("Standard acceleration of gravity", 9.80665e-14,
                   pc.accelerationFromGravity(), pc.g(), ut);
  compare_constant("Faraday constant", 8.667420510246794e+16,
                   pc.faradayConstant(), pc.F(), ut);
  compare_constant("Permeability of vacuum", 1.25663706212e-15,
                   pc.permeabilityOfVacuum(), pc.mu0(), ut);
  compare_constant("Permittivity of vacuum", 0.07953848254861551,
                   pc.permittivityOfVacuum(), pc.epsi0(), ut);
  compare_constant("Classical electron radius", 2.8179403262e-13,
                   pc.classicalElectronRadius(), pc.re(), ut);
  compare_constant("Electron mass", 9.1093837015e-28, pc.electronMass(),
                   pc.me(), ut);
  compare_constant("Proton mass", 1.67262192369e-24, pc.protonMass(), pc.mp(),
                   ut);
  compare_constant("Electron volt", 1.602176634e-28, pc.electronVolt(), pc.eV(),
                   ut);
}

} // namespace rtt_units_test

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    rtt_units_test::test_si(ut);
    rtt_units_test::test_cgs(ut);
    rtt_units_test::test_cgsh(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// End of test/tCompton.cc
//----------------------------------------------------------------------------//
