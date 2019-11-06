//----------------------------------*-C++-*----------------------------------//
/*!                                                                              
 * \file   cdi_CPEloss/Tabular_CP_Eloss.cc                                     
 * \author Ben R. Ryan                                                        
 * \date   2019 Nov 4                                             
 * \brief  Tabular_CP_Eloss member definitions.                                 
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.                 
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Tabular_CP_Eloss.hh"

using std::stod;
using std::stoi;

namespace rtt_cdi_cpeloss {

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic tabular eloss model. 
 *
 * This constructor builds an eloss model defined by the
 * rtt_cdi_cpeloss::Tabular_Eloss_Model derived class argument.
 *
 * The path to an eloss datafile is passed to the constructor,
 * which opens and parses the file. The file format is the usual
 * LANL format for stopping powers. 
 *
 * \param filename_in path to eloss file
 * \param target_zaid_in target particle zaid
 * \param projectile_zaid_in transporting particle zaid
 */
Tabular_CP_Eloss::Tabular_CP_Eloss(std::string filename_in,
                                   rtt_cdi::CParticle &target_in,
                                   rtt_cdi::CParticle &projectile_in)
    : rtt_cdi::CPEloss(target_in, projectile_in) {
    //: filename(filename_in) rtt_cdi::CPEloss(target_in, projectile_in) {
  filename = filename_in;
        
  model_type = rtt_cdi::CPModelType::TABULAR_ETYPE;

  file.open(filename);
  if (!file.is_open()) {
    std::cout << "Eloss file " << filename << " could not be opened!"
              << std::endl;
    exit(-1);
  }

  std::string line;
  std::vector<std::string> line_entries;
  int nlines;
  int max_entries =
      6; // This is a statement about the file format, maximum of six entries for row.

  read_line(); // ZAID
  int32_t projectile_zaid_file = stoi(line_entries[0]);
  Require(projectile.get_zaid() == projectile_zaid_file);

  read_line(); // Z, A, mass

  line_entries = read_line(); // Number of bins for energy, density, temperature
  n_energy = stoi(line_entries[0]);
  n_density = stoi(line_entries[1]);
  n_temperature = stoi(line_entries[2]);

  line_entries =
      read_line(); // Bin spacing for energy, density, temperature (log)
  d_log_energy = 1. / stod(line_entries[0]);
  d_log_density = 1. / stod(line_entries[1]);
  d_log_temperature = 1. / stod(line_entries[2]);

  // Get first energy support point
  nlines = std::ceil((double)n_energy / max_entries);
  line_entries = read_line();
  min_log_energy = stod(line_entries[0]);
  for (int n = 0; n < nlines - 1; n++) {
    read_line();
  }

  // Get first density support point
  nlines = std::ceil((double)n_density / max_entries);
  line_entries = read_line();
  min_log_density = stod(line_entries[0]);
  for (int n = 0; n < nlines - 1; n++) {
    read_line();
  }

  // Get first temperature support point
  nlines = std::ceil((double)n_temperature / max_entries);
  line_entries = read_line();
  min_log_temperature = stod(line_entries[0]);
  for (int n = 0; n < nlines - 1; n++) {
    read_line();
  }

  std::vector<double> stopping_data_1d(n_energy * n_density * n_temperature);

  bool target_found = false;
  nlines = std::ceil(
      ((double)n_energy * n_density * n_temperature) /
      max_entries); // The number of lines taken up by stopping power data for one target
  if (target.get_zaid() == -1) {
    // Target is free electrons
    target_found = true;
    //nlines = std::ceil(((double)n_energy*n_density*n_temperature)/max_entries);
    int nentry = 0;
    for (int n = 0; n < nlines; n++) {
      line_entries = read_line();
      for (std::string entry : line_entries) {
        stopping_data_1d[nentry] = stod(entry);
        nentry++;
      }
    }
  } else {
    // Skip electrons
    //nlines = std::ceil(((double)n_energy*n_density*n_temperature)/max_entries);
    for (int n = 0; n < nlines; n++) {
      line_entries = read_line();
    }

    // Find ion target, if it exists
    int n_target_ions = stoi(read_line()[0]); // Number of target ions in file
    for (int n_target_ion = 0; n_target_ion < n_target_ions; n_target_ion++) {
      int zaid_target_ion = stoi(read_line()[0]); // ZAID
      read_line();                                // Z, A, mass
      if (zaid_target_ion == target.get_zaid()) {
        // This is the requested target ion
        target_found = true;
        int nentry = 0;
        for (int n = 0; n < nlines; n++) {
          line_entries = read_line();
          for (std::string entry : line_entries) {
            stopping_data_1d[nentry] = stod(entry);
            nentry++;
          }
        }
        break;
      } else {
        // This is not the requested target ion
        for (int n = 0; n < nlines; n++) {
          read_line();
        }
      }
    }
  }
  file.close();

  if (!target_found) {
    std::cout << "Target ZAID " << target.get_zaid() << " not found in DEDX fle "
              << filename << "!" << std::endl;
    exit(-1);
  }

  rtt_dsxx::DracoArray<double> stopping_data(n_energy, n_density,
                                             n_temperature);
  for (int ne = 0; ne < n_energy; ne++) {
    for (int nd = 0; nd < n_density; nd++) {
      for (int nt = 0; nt < n_temperature; nt++) {
        stopping_data(nt, nd, ne) =
            stopping_data_1d[ne + n_energy * (nd + n_density * nt)];
      }
    }
  }

  // Convert units on table to match those of getEloss:
  //   energy:      MeV -> cm/shk (using target particle mass)
  double energy = exp(min_log_energy)*(1.e6*1.6021772e-12);
  min_log_energy = log(sqrt(2.*energy/target.get_mass())*1.e-8);
  d_log_energy = d_log_energy / 2.;
  //   density:     cm^-3 -> g cm^-3
  min_log_density = log(exp(min_log_density)*target.get_mass());
  //   temperature: keV -> keV
  // Note that d log x = dx / x is not affected by linear unit conversions

  // Initialize table bounds
  min_energy = exp(min_log_energy);
  max_energy = exp(min_log_energy + d_log_energy*n_energy);
  min_density = exp(min_log_density);
  max_density = exp(min_log_density + d_log_density*n_density);
  min_temperature = exp(min_log_temperature);
  max_temperature = exp(min_log_temperature + d_log_temperature*n_temperature);
}

/*!
 * \brief Read a line from an eloss datafile and return as a vector of strings.
 *
 * Convenience function, especially when reading header data which is not
 * uniform in the number of entries per line, or the types of those entries.
 *
 * \return entries the resulting vector of entries in the datafile line.
 *
 */
std::vector<std::string> Tabular_CP_Eloss::read_line() {
  std::string line;
  std::getline(file, line);
  std::vector<std::string> entries;
  std::istringstream iss(line);
  for (std::string s; iss >> s;) {
    entries.push_back(s);
  }
  return entries;
}

double Tabular_CP_Eloss::getEloss(const double temperature, const double density, const double partSpeed) const {
  Require(temperature >= 0.);
  Require(density >= 0.);
  Require(partSpeed >= 0.);

  if (temperature < min_temperature || temperature > max_temperature ||
      density < min_density || density > max_density ||
      partSpeed < min_energy || partSpeed > max_energy) {
    // Outside of the table
    return 0.;
  }

  return 0.;
  //return rtt_dsxx::linear_interpolate_3(x0, x1, y0, y1, z0, z1, f000, f100, f001, f101, f010, f110, f011, f111, x, y, z);
}

} // namespace rtt_cdi_cpeloss
