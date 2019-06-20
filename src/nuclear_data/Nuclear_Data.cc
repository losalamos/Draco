//----------------------------------*-C++-*--------------------------------//
/*!
 * \file   nuclear_data/Nuclear_Data.cc
 * \author B. R. Ryan
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for nuclear_data/Nuclear_Data class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Nuclear_Data.hh"

namespace rtt_nuclear_data {

//----------------------------------------------------------------------------//
/*!
 * \brief Parses the cell_data block data from the mesh file via calls to
 *        private member functions. See Appendix F of the MCNP User Manual for
 *        details on the ACE file format.
 * \param _filepath Path to ACE file.
 */
Nuclear_Data::Nuclear_Data(const string &_filepath) {
  filepath = _filepath;
  string filename = filepath.substr(filepath.find_last_of("/") + 1);
  Insist(filename.length() <= 10, "Filename has nonstandard length");

  const char *test = filepath.c_str();
  ifstream ACEfile(test, std::ios::in);
  if (!ACEfile) {
    std::ostringstream buffer;
    buffer << "File " << filename << " could not be opened\n";
    throw std::invalid_argument(buffer.str());
  }

  string line;
  std::getline(ACEfile, line);
  zaid = line.substr(0, 10);
  atomic_weight = stof(line.substr(10, 12));
  temperature = stof(line.substr(22, 12));
  date = line.substr(24, 10);

  std::getline(ACEfile, line);
  comment = line.substr(0, 70);
  mat_identifier = line.substr(70, 10);

  string reaction_name = zaid.substr(9, 1);
  if (!reaction_name.compare("c") || !reaction_name.compare("C")) {
    reaction = Reaction::CONTINUOUS_ENERGY_NEUTRON;
  } else if (!reaction_name.compare("d") || !reaction_name.compare("D")) {
    reaction = Reaction::DISCRETE_REACTION_NEUTRON;
  } else if (!reaction_name.compare("y") || !reaction_name.compare("Y")) {
    reaction = Reaction::DOSIMETRY;
  } else if (!reaction_name.compare("t") || !reaction_name.compare("T")) {
    reaction = Reaction::THERMAL_S;
  } else if (!reaction_name.compare("p") || !reaction_name.compare("P")) {
    reaction = Reaction::CONTINUOUS_ENERGY_PHOTOATOMIC;
  } else if (!reaction_name.compare("m") || !reaction_name.compare("M")) {
    reaction = Reaction::NEUTRON_MULTIGROUP;
  } else if (!reaction_name.compare("g") || !reaction_name.compare("G")) {
    reaction = Reaction::PHOTOATOMIC_MULTIGROUP;
  } else if (!reaction_name.compare("e") || !reaction_name.compare("E")) {
    reaction = Reaction::CONTINUOUS_ENERGY_ELECTRON;
  } else if (!reaction_name.compare("u") || !reaction_name.compare("U")) {
    reaction = Reaction::CONTINUOUS_ENERGY_PHOTONUCLEAR;
  }

  atomic_number = stoi(zaid.substr(0,3));
  mass_number = stoi(zaid.substr(3,3)); // 000 for P, M, G, E
  evaluation_identifier = stoi(zaid.substr(7, 2));

  for (int n = 0; n < 4; n++) {
    std::getline(ACEfile, line);
    for (int m = 0; m < 4; m++) {
      IZ.push_back(stoi(line.substr(0+18*m,7)));
      AW.push_back(stof(line.substr(7+18*m,11)));
    }
  }

  string NXS_s, JXS_s;
  for (int n = 0; n < 2; n++) {
    std::getline(ACEfile, line);
    NXS_s.append(line);
  }
  for (int n = 0; n < 4; n++) {
    std::getline(ACEfile, line);
    JXS_s.append(line);
  }
  for (int n = 0; n < 16; n++) {
    NXS.push_back(stoi(NXS_s.substr(0+9*n,9)));
  }
  std::cout << NXS_s << std::endl;
  std::cout << JXS_s << std::endl;
  for (int n = 0; n < 16; n++) {
    std::cout << NXS[n] << std::endl;
  }

  // Interpret file depending on reaction type
  if (reaction == Reaction::CONTINUOUS_ENERGY_NEUTRON || 
      reaction == Reaction::DISCRETE_REACTION_NEUTRON) {
    // NXS
    data_length = NXS[0];
    ZA = NXS[1];
  } else if (reaction == Reaction::DOSIMETRY) {
    return;
  } else if (reaction == Reaction::THERMAL_S) {
    return;
  } else if (reaction == Reaction::CONTINUOUS_ENERGY_PHOTOATOMIC) {
    data_length = NXS[0];
  } else if (reaction == Reaction::NEUTRON_MULTIGROUP) {
    return;
  } else if (reaction == Reaction::PHOTOATOMIC_MULTIGROUP) {
    return;
  } else if (reaction == Reaction::CONTINUOUS_ENERGY_ELECTRON) {
    return;
  } else if (reaction == Reaction::CONTINUOUS_ENERGY_PHOTONUCLEAR) {
    return;
  }
}

void Nuclear_Data::report_contents() {
  std::cout << "File:                  " << filepath << std::endl;
  std::cout << "ZAID:                  " << zaid << std::endl;
  std::cout << "Atomic Weight:         " << atomic_weight << std::endl;
  std::cout << "Temperature:           " << temperature << std::endl;
  std::cout << "Date:                  " << date << std::endl;
  std::cout << "Comment:               " << comment << std::endl;
  std::cout << "MAT Identifier:        " << mat_identifier << std::endl;
  std::cout << "Atomic Number:         " << atomic_number << std::endl;
  std::cout << "Mass Number:           " << mass_number << std::endl;
  std::cout << "Evaluation Identifier: " << evaluation_identifier << std::endl;
  std::cout << "Reaction:              " << Reaction_Names[reaction] << std::endl;
  std::cout << "IZ, AW pairs:          " << std::endl;
  for (int n = 0; n < 16; n++) {
    std::cout << "  " << IZ[n] << " " << AW[n] << std::endl;
  }
  std::cout << "NXS" << std::endl;
  std::cout << "data_length: " << data_length << std::endl;
}
} // end namespace rtt_nuclear_data
