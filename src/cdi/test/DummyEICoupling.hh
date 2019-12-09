//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyEICoupling.hh
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  DummyEICoupling class header file (derived from ../EICoupling)
 * \note   Copyright (C) 2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyEICoupling_hh__
#define __cdi_DummyEICoupling_hh__

#include "cdi/EICoupling.hh"

namespace rtt_cdi_test {

//========================================================================
/*!
 * \class DummyEICoupling
 *
 * \brief This is an equation of state class that derives its interface from
 *        cdi/EICoupling and is used for testing purposes only.
 *
 * This electron-ion coupling class always contains the same data (set by the
 * default constructor).  The data table has the following properties:
 *
 *    dummy_ei_coupling = etemperature + 10*itemperature + 100*density +
 *    1000*w_e + 10000*w_i 
 *
 * \sa cdi/test/tEICoupling.cc
 * \sa cdi/test/tCDI.cc
 */
//========================================================================

class DLL_PUBLIC_cdi_test DummyEICoupling : public rtt_cdi::EICoupling {
public:
  // -------------------------- //
  // Constructors & Destructors //
  // -------------------------- //

  /*!
   * \brief Constructor for DummyEICoupling object.
   *
   * The constructor assigns fixed values for all of the member data.  Every
   * instance of this object has the same member data.
   */
  DummyEICoupling();

  /*!
   * \brief Default DummyEICoupling() destructor.
   *
   * This is required to correctly release memory when a DummyEICoupling object is
   * destroyed.
   */
  ~DummyEICoupling();

  // --------- //
  // Accessors //
  // --------- //

  /*!
   * \brief EICoupling accessor that returns a single electron-ion coupling
   * given an electron and ion temperature, the material density, and the
   * electron and ion screening coeffiecients..
   *
   *    dummy_ei_coupling = etemperature + 10*itemperature + 100*density +
   *    1000*w_e + 10000*w_i
   *
   * \param[in] etemperature The electron temperature value for which an opacity
   *        value is being requested (Kelvin).
   * \param[in] itemperature The ion temperature value for which an opacity
   *        value is being requested (Kelvin).
   * \param[in] density The density value for which an opacity value is being
   *        requested (g/cm^3).
   * \param[in] w_e the electron screening coeffiecent [1/s]
   * \param[in] w_i the ion screening coeffiecent [1/s]
   * \return An electron-ion coupling coeffient [kJ/g/K/s].
   */
  double getElectronIonCoupling(const double etemperature,
                                const double itemperature, const double density,
                                const double w_e, const double w_i) const;

  /*!
   * \brief EICoupling accessor that returns a vector of electron-ion coupling
   * given an electron and ion temperature, the material density, and the
   * electron and ion screening coeffiecients..
   *
   *    dummy_ei_coupling = etemperature + 10*itemperature + 100*density +
   *    1000*w_e + 10000*w_i
   *
   * \param[in] etemperature The electron temperature vector for which an opacity
   *        value is being requested (Kelvin).
   * \param[in] itemperature The ion temperature vector for which an opacity
   *        value is being requested (Kelvin).
   * \param[in] density The density vector for which an opacity value is being
   *        requested (g/cm^3).
   * \param[in] w_e the electron screening coeffiecent vector [1/s]
   * \param[in] w_i the ion screening coeffiecent vector [1/s]
   * \return An electron-ion coupling coeffient vector [kJ/g/K/s].
   */
  std::vector<double>
  getElectronIonCoupling(const std::vector<double> &etemperature,
                         const std::vector<double> &itemperature,
                         const std::vector<double> &density,
                         const std::vector<double> &w_e,
                         const std::vector<double> &w_i) const;

  // Dummy pack function.
  std::vector<char> pack() const { return std::vector<char>(); }

}; // end of class DummyEICoupling

} // end namespace rtt_cdi_test

#endif // __cdi_DummyEICoupling_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/DummyEICoupling.hh
//---------------------------------------------------------------------------//
