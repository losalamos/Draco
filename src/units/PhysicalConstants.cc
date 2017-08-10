//----------------------------------*-C++-*----------------------------------//
/*! \file   PhysicalConstants.cc
 *  \author Kelly Thompson
 *  \brief  Provide a single place where physical constants (pi, speed of
 *          light, etc) are defined for the local UnitSystem.
 *  \date   Mon Nov 10 09:24:55 2003
 *  \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *          All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "PhysicalConstants.hh"

#include "PhysicalConstantsSI.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iomanip>
#include <iostream>

namespace rtt_units {

/*!
 * \brief Default constructor provides physical constants with SI units (kg,
 *        m, seconds, degree K, amp, radian, mole).
 * \return A PhysicalConstants object.
 */
PhysicalConstants::PhysicalConstants()
    : d_avogadro(AVOGADRO), d_planck(planckSI), d_gasConstant(gasConstantSI),
      d_boltzmann(boltzmannSI), d_electronCharge(electronChargeSI),
      d_cLight(cLightSI), d_stefanBoltzmann(stefanBoltzmannSI),
      d_gravitationalConstant(gravitationalConstantSI),
      d_accelerationFromGravity(accelerationFromGravitySI),
      d_faradayConstant(faradayConstantSI),
      d_permeabilityOfVacuum(permeabilityOfVacuumSI),
      d_permittivityOfFreeSpace(permittivityOfFreeSpaceSI),
      d_classicalElectronRadius(classicalElectronRadiusSI),
      d_electronMass(electronMassSI), d_protonMass(protonMassSI) {
  // empty
}

/*!
 * \brief Constuctor that creates a set of PhysicalConstants using the
 *        provided rtt_units::UnitSystem.
 * \param u A complete UnitSystem.  PhysicalConstants will be formed and
 *        returned using these units (CGS, SI, etc.)
 * \return A PhysicalConstants object.
 */
PhysicalConstants::PhysicalConstants(UnitSystem const &u)
    : d_avogadro(AVOGADRO * u.Q()), d_planck(planckSI * u.e() * u.t()),
      d_gasConstant(gasConstantSI * u.e() / u.T()),
      d_boltzmann(boltzmannSI * u.e() / u.T()),
      d_electronCharge(electronChargeSI), d_cLight(cLightSI * u.v()),
      d_stefanBoltzmann(stefanBoltzmannSI * u.p() / std::pow(u.L(), 2) /
                        std::pow(u.T(), 4)),
      d_gravitationalConstant(gravitationalConstantSI * u.f() *
                              std::pow(u.L(), 2) / std::pow(u.M(), 2)),
      d_accelerationFromGravity(accelerationFromGravitySI * u.a()),
      d_faradayConstant(AVOGADRO * d_electronCharge),
      d_permeabilityOfVacuum(permeabilityOfVacuumSI / u.L()),
      d_permittivityOfFreeSpace(1.0 / d_permeabilityOfVacuum /
                                std::pow(d_cLight, 2)),
      d_classicalElectronRadius(classicalElectronRadiusSI * u.L()),
      d_electronMass(electronMassSI * u.M()),
      d_protonMass(protonMassSI * u.M()) {
  // empty

  // std::cout << std::setprecision(16) << std::scientific
  //     << "\nu.e() = " << u.e()
  //     << "\nu.t() = " << u.t()
  //     << "\nu.T() = " << u.T()
  //     << "\nu.v() = " << u.v()
  //     << "\nu.p() = " << u.p()
  //     << "\nu.L() = " << u.L()
  //     << "\nu.M() = " << u.M()
  //     << "\nu.a() = " << u.a()
  //     << std::endl;
}

} // end namespace rtt_units

//---------------------------------------------------------------------------//
// end of units/PhysicalConstants.cc
//---------------------------------------------------------------------------//
