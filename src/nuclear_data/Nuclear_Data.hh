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
const char* Reaction_Names[] = { "CONTINUOUS_ENERGY_NEUTRON", 
  "DISCRETE_REACTION_NEUTRON", "DOSIMETRY", "THERMAL_S", 
  "Continuous-energy photoatomic", "NEUTRON_MULTIGROUP", 
  "PHOTOATOMIC_MULTIGROUP", "CONTINUOUS_ENERGY_ELECTRON", 
                      "CONTINUOUS_ENERGY_PHOTONUCLEAR" };

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

  // NXS array
  // Not all values are used for each reaction type; see Table F.1 of MCNP User 
  // Manual.
  std::vector<int> NXS;
  int data_length;       // Length of primary data block
  int ZA;                // 1000*Z + A
  int IDPNI;             // Inelastic scattering mode
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

  // JXS array
  // Not all values are used for each reaction type; see Table F.2 of MCNP User
  // Manual.
  std::vector<int> JXS;
  int ESZ;   // Location of energy table
  int LONE;  // Location of first word of table
  int ITIE;  // Location of inelastic energy table
  int ESZG;  // Location of energy table
  int NU;    // Location of fission nu data
  int ITIX;  // Location of inelastic cross sections
  int JINC;  // Location of incoherent form factors
  int MTR;   // Location of MT array
  int ITXE;  // Location of inelastic energy/angle distributions
  int JCOH;  // Location of coherent form factors
  int LQR;   // Location of Q-value array
  int ITCE;  // Location of elastic energy table
  int JFLO;  // Location of fluorescence data
  int TYR;   // Location of reaction type array
  int ITCX;  // Location of elastic cross sections
  int LHNM;  // Location of heating numbers
  int LSIG;  // Location of cross-section locators
  int ITCA;  // Location of elastic angular distributions
  int LNEPS; // Location of the number of electrons per shell
  int SIG;   // Location of cross sections
  int SIGD;  // Location of cross sections
  int LBEPS; // Location of binding energy per shell
  int LAND;  // Location of table of angular distribution locators
  int LPIPS; // Location of probability of interaction per shell
  int AND;   // Location of angular distributions
  int LSWD;  // Location of array of offsets to shellwise data
  int LDLW;  // Location of table of energy distribution locators
  int SWD;   // Location of shellwise data in PDF and CDF form
  int DLW;   // Location of energy distributions
  int GPD;   // Location of photon production data
  int MTRP;  // Location of photon production MT array
  int LSIGP; // Location of table of photon production cross-section locators
  int SIGP;  // Location of photon production cross sections
  int LANDP; // Location of table of photon production angular distribution
             // locators
  int ANDP;  // Location of photon production angular distributions
  int LDLWP; // Location of table of photon production energy distribution
             // locators
  int DLWP;  // Location of photon production energy distributions
  int YP;    // Location of table of yield multipliers
  int FIS;   // Location of total fission cross section
  int END;   // Location of last word of this table
  int LUNR;  // Location of probability tables
  int DNU;   // Location of delayed nubar data
  int BDD;   // Location of basic delayed data (lambdas, probabilities)
  int DNEDL; // Location of table of energy distribution locators
  int DNED;  // Location of energy distributions
};
} // end namespace rtt_nuclear_data

#endif // __nuclear_data_Nuclear_Data_hh__

//---------------------------------------------------------------------------//
// end of nuclear_data/Nuclear_Data.hh
//---------------------------------------------------------------------------//
