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
Nuclear_Data::Nuclear_Data(const string &filename) {
  Insist(filename.length() <= 10, "Filename has nonstandard length");

  const char *test = filename.c_str();
  ifstream ACEfile(test, std::ios::in);
  if (!ACEfile) {
    std::ostringstream buffer;
    buffer << "File " << filename << " could not be opened\n";
    throw std::invalid_argument(buffer.str());
  }

  string reaction_name = filename.substr(filename.find(".") + 1).substr(2);
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
  std::cout << reaction << std::endl;
}
} // end namespace rtt_nuclear_data
