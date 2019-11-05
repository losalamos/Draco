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
#include "ds++/DracoArray.hh"

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
                                   int32_t target_zaid_in,
                                   int32_t projectile_zaid_in)
    : filename(filename_in), target_zaid(target_zaid_in),
      projectile_zaid(projectile_zaid_in) {
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
  Require(projectile_zaid == projectile_zaid_file);

  read_line(); // Z, A, mass
  
  line_entries = read_line(); // Number of bins for energy, density, temperature
  n_energy = stoi(line_entries[0]);
  n_density = stoi(line_entries[1]);
  n_temperature = stoi(line_entries[2]);
  
  line_entries = read_line(); // Bin spacing for energy, density, temperature (log)
  d_log_energy = 1./stod(line_entries[0]);
  d_log_density = 1./stod(line_entries[1]);
  d_log_temperature = 1./stod(line_entries[2]);
  
  // Get first energy support point
  nlines = std::ceil((double)n_energy/max_entries);
  line_entries = read_line();
  min_log_energy = stod(line_entries[0]);
  for (int n = 0; n < nlines - 1; n++) {
    read_line();
  }

  // Get first density support point
  nlines = std::ceil((double)n_density/max_entries);
  line_entries = read_line();
  min_log_density = stod(line_entries[0]);
  for (int n = 0; n < nlines - 1; n++) {
    read_line();
  }
  
  // Get first temperature support point
  nlines = std::ceil((double)n_temperature/max_entries);
  line_entries = read_line();
  min_log_temperature = stod(line_entries[0]);
  for (int n = 0; n < nlines - 1; n++) {
    read_line();
  }

  rtt_dsxx::DracoArray<double> stopping_data(n_energy, n_density, n_temperature);
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
  getline(file, line);
  std::vector<std::string> entries;
  std::istringstream iss(line);
  for (std::string s; iss >> s;) {
    entries.push_back(s);
  }
  return entries;
}

} // namespace rtt_cdi_cpeloss
