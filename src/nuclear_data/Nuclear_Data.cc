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
 *        private member functions.
 * \param meshfile Mesh file name.
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

  //zaid = string(10 - filename.length(), '0').append(filename);
  zaid = line.substr(0, 10);

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

  // Interpret file depending on reaction type
  if (reaction == Reaction::CONTINUOUS_ENERGY_NEUTRON) {
  } else if (reaction == Reaction::CONTINUOUS_ENERGY_PHOTOATOMIC) {
    data_length = 0;
  }

  string x;
  std::getline(ACEfile, x);
  std::cout << x << std::endl;

}

void Nuclear_Data::report_contents() {
  std::cout << "File:                  " << filepath << std::endl;
  std::cout << "ZAID:                  " << zaid << std::endl;
  std::cout << "Atomic Number:         " << atomic_number << std::endl;
  std::cout << "Mass Number:           " << mass_number << std::endl;
  std::cout << "Evaluation Identifier: " << evaluation_identifier << std::endl;
  std::cout << "Reaction:              " << reaction << std::endl;
}
} // end namespace rtt_nuclear_data
