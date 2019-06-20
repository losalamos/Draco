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
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;

public:
  Nuclear_Data(const string &filename);
  ~Nuclear_Data() {}

  void report_contents();

private:
  string filepath;
  string zaid;
  double atomic_weight;
  double temperature;
  string date;
  string comment;
  string mat_identifier;

  int atomic_number;
  int mass_number;
  string thermal_abbreviation;
  int evaluation_identifier;
  int reaction;

  // Z/AW pairs, now only used for thermal tables to indicate applicable 
  // isotopes.
  std::vector<int> IZ;
  std::vector<double> AW;

  // NXS Array
  // Not all values are used for each reaction type; see Table F.1 of MCNP User 
  // Manual.
  int data_length;       // Length of primary data block
  int ZA;                // 1000*Z + A
  int IDPNI;             // Inelastic scatteirng mode
  int Z;                 // Z
  int NES;               // Number of energies
  int NIL;               // Inelastic dimensioning parameter
  int NTR;               // Number of reactions excluding elastic
  int NIEB;              // Number of inelastic exiting energies
  int NFLO;              // Length of the fluorescence data divided by 4
  int NR;                // Number of reactions having secondary neutrons 
                         // excluding elastic
  int IDPNC;             // Elastic scattering mode
  int NSH;               // Number of electron shells
  int NTRP;              // Number of photon production reactions
  int NCL;               // Elastic dimensioning parameter
  int IFENG;             // Secondary energy mode
  int NPCR;              // Number of delayed neutron precursor families
  int NT;                // Number of PIKMT reactions
  int photon_production; //  0 = normal photon production
                         // -1 = do not produce photons

  int data_length;
};
} // end namespace rtt_nuclear_data

#endif // __nuclear_data_Nuclear_Data_hh__

//---------------------------------------------------------------------------//
// end of nuclear_data/Nuclear_Data.hh
//---------------------------------------------------------------------------//
