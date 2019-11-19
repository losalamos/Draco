//----------------------------------*-C++-*----------------------------------//
/*!                                                                              
 * \file   cdi_CPEloss/Tabular_CP_Eloss.hh                                     
 * \author Ben R. Ryan                                                        
 * \date   2019 Nov 4                                              
 * \brief  Tabular_CP_Eloss class definition.                                   
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.                 
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_CPEloss_Tabular_CP_Eloss_hh__
#define __cdi_CPEloss_Tabular_CP_Eloss_hh__

#include "cdi/CPCommon.hh"
#include "cdi/CPEloss.hh"
#include "ds++/Assert.hh"
#include "ds++/DracoArray.hh"
#include "ds++/DracoMath.hh"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace rtt_cdi_cpeloss {
//===========================================================================//
/*!                                                                              
 * \class Tabular_CP_Eloss                                                      
 *                                                                               
 * \brief Derived rtt_cdi::CPEloss class for tabular eloss.                     
 *        This class implements the interface found in cdi/CPEloss.hh for        
 *        the case where CP energy loss data is in tabular form, stored 
 *        in a file.                
 *                                                                               
 */
//===========================================================================//

class Tabular_CP_Eloss : public rtt_cdi::CPEloss {
public:
typedef std::vector<double> sf_double;

private:
  std::string filename;
  std::ifstream file;

  uint32_t n_energy;           // Number of gridpoints in projectile energy
  uint32_t n_density;          // Number of gridpoints in target density
  uint32_t n_temperature;      // Number of gridpoints in target temperature
  double d_log_energy;        // Log spacing of projectile energy gridpoints
  double d_log_density;       // Log spacing of target density gridpoints
  double d_log_temperature;   // Log spacing of target temperature gridpoints
  double min_log_energy;      // Log of minimum projectile energy
  double min_log_density;     // Log of minimum target density
  double min_log_temperature; // Log of minimum target temperature
  double min_energy;          // Minimum target energy
  double max_energy;          // Maximum target energy
  double min_density;         // Minimum target density
  double max_density;         // Maximum target density
  double min_temperature;     // Minimum target temperature
  double max_temperature;     // Maximum target temperature
  sf_double energies;         // Vector of energy gridpoints
  sf_double densities;        // Vector of density gridpoints
  sf_double temperatures;     // Vector of temperature gridpoints
  // Note that after unit conversions, *_energy is really *_speed

  // Storage for tabulated data
  rtt_dsxx::DracoArray<double> stopping_data;

  // Utility for reading a line of an eloss file and as a vector of strings
  std::vector<std::string> read_line();

public:
  // Constructor
  Tabular_CP_Eloss(
      std::string filename_in, rtt_cdi::CParticle target_in,
      rtt_cdi::CParticle
          projectile_in);

  // >>> ACCESSORS

  double getEloss(const double temperature, const double density,
                  const double v0) const;

  bool data_in_tabular_form() const { return true; }

  inline std::string getDataFilename() const { return filename; }

  sf_double getTemperatureGrid() const { return temperatures; }

  sf_double getDensityGrid() const { return densities; }

  sf_double getEnergyGrid() const { return energies; }

  size_t getNumTemperatures() const { return n_temperature; }

  size_t getNumDensities() const { return n_density; }

  size_t getNumEnergies() const { return n_energy; }
};

} // namespace rtt_cdi_cpeloss

#endif
