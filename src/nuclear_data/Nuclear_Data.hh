//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   nuclear_data/Nuclear_Data.hh
 * \author B. R. Ryan
 * \date   Wed Jun 19 2019
 * \brief  Header file for Nuclear_Data class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __nuclear_data_Nuclear_Data_hh__
#define __nuclear_data_Nuclear_Data_hh__

#include "ds++/Assert.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace rtt_nuclear_data {
enum Reaction { CONTINUOUS_ENERGY_NEUTRON, DISCRETE_REACTION_NEUTRON,
                DOSIMETRY, THERMAL_S, CONTINUOUS_ENERGY_PHOTOATOMIC,
                NEUTRON_MULTIGROUP, PHOTOATOMIC_MULTIGROUP,
                CONTINUOUS_ENERGY_ELECTRON, CONTINUOUS_ENERGY_PHOTONUCLEAR };

//===========================================================================//
/*!
 * \class Nuclear_Data
 *
 * \brief A class for storing the data from one ACE nuclear data file.
 *
 *\sa The Nuclear_Data class contains data members for storing all the data in
 *    an ACE nuclear data file. The name of an ACE file is passed to the class
 *    constructor, which then reads that file and populates the class's data
 *    members. Public functions are provided for accessing this data.
 */
//===========================================================================//
class Nuclear_Data {
//class Nuclear_Data {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;

public:
  Nuclear_Data(const string &filename);
  ~Nuclear_Data() {}

private:
  int atomic_number;
  int mass_number;
  string thermal_abbreviation;
  int evaluation_identifier;
  int reaction;
};
} // end namespace rtt_nuclear_data

#endif // __nuclear_data_Nuclear_Data_hh__

//---------------------------------------------------------------------------//
// end of nuclear_data/Nuclear_Data.hh
//---------------------------------------------------------------------------//
