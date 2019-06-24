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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace rtt_nuclear_data {
enum Reaction {
  CONTINUOUS_ENERGY_NEUTRON,
  DISCRETE_REACTION_NEUTRON,
  DOSIMETRY,
  THERMAL_S,
  CONTINUOUS_ENERGY_PHOTOATOMIC,
  NEUTRON_MULTIGROUP,
  PHOTOATOMIC_MULTIGROUP,
  CONTINUOUS_ENERGY_ELECTRON,
  CONTINUOUS_ENERGY_PHOTONUCLEAR
};
const char *Reaction_Name[] = {"Continuous-energy neutron",
                               "Discrete-reaction neutron",
                               "Dosimetry",
                               "Thermal S(alpha, beta)",
                               "Continuous-energy photoatomic",
                               "Neutron multigroup",
                               "Photoatomic multigroup",
                               "Continuous-energy electron",
                               "Continuous-energy photonuclear"};

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

private:
  struct Datatable {
    int start_index = -1;
    int length;
    std::vector<double> data;

    void set_indices(int _start_index, int _length) {
      start_index = _start_index;
      length = _length;
    }
  };

public:
  Nuclear_Data(const string &_filepath);
  ~Nuclear_Data() {}

  void report_contents();

  const std::vector<double> &get_ESZ() { return ESZ.data; };
  const std::vector<double> &get_LONE() { return LONE.data; };
  const std::vector<double> &get_ITIE() { return ITIE.data; };
  const std::vector<double> &get_ESZG() { return ESZG.data; };
  const std::vector<double> &get_NU() { return NU.data; };
  const std::vector<double> &get_ITIX() { return ITIX.data; };
  const std::vector<double> &get_JINC() { return JINC.data; };
  const std::vector<double> &get_MTR() { return MTR.data; };
  const std::vector<double> &get_ITXE() { return ITXE.data; };
  const std::vector<double> &get_JCOH() { return JCOH.data; };
  const std::vector<double> &get_LQR() { return LQR.data; };
  const std::vector<double> &get_ITCE() { return ITCE.data; };
  const std::vector<double> &get_JFLO() { return JFLO.data; };
  const std::vector<double> &get_TYR() { return TYR.data; };
  const std::vector<double> &get_ITCX() { return ITCX.data; };
  const std::vector<double> &get_LHNM() { return LHNM.data; };
  const std::vector<double> &get_LSIG() { return LSIG.data; };
  const std::vector<double> &get_ITCA() { return ITCA.data; };
  const std::vector<double> &get_LNEPS() { return LNEPS.data; };
  const std::vector<double> &get_SIG() { return SIG.data; };
  const std::vector<double> &get_SIGD() { return SIGD.data; };
  const std::vector<double> &get_LBEPS() { return LBEPS.data; };
  const std::vector<double> &get_LAND() { return LAND.data; };
  const std::vector<double> &get_LPIPS() { return LPIPS.data; };
  const std::vector<double> &get_AND() { return AND.data; };
  const std::vector<double> &get_LSWD() { return LSWD.data; };
  const std::vector<double> &get_LDLW() { return LDLW.data; };
  const std::vector<double> &get_SWD() { return SWD.data; };
  const std::vector<double> &get_DLW() { return DLW.data; };
  const std::vector<double> &get_GPD() { return GPD.data; };
  const std::vector<double> &get_MTRP() { return MTRP.data; };
  const std::vector<double> &get_LSIGP() { return LSIGP.data; };
  const std::vector<double> &get_SIGP() { return SIGP.data; };
  const std::vector<double> &get_LANDP() { return LANDP.data; };
  const std::vector<double> &get_ANDP() { return ANDP.data; };
  const std::vector<double> &get_LDLWP() { return LDLWP.data; };
  const std::vector<double> &get_DLWP() { return DLWP.data; };
  const std::vector<double> &get_YP() { return YP.data; };
  const std::vector<double> &get_FIS() { return FIS.data; };
  const std::vector<double> &get_END() { return END.data; };
  const std::vector<double> &get_LUNR() { return LUNR.data; };
  const std::vector<double> &get_DNU() { return DNU.data; };
  const std::vector<double> &get_BDD() { return BDD.data; };
  const std::vector<double> &get_DNEDL() { return DNEDL.data; };
  const std::vector<double> &get_DNED() { return DNED.data; };

private:
  void load_datatable(Datatable &dataTable, const std::string data_s);

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

  // JXS array, with associated data
  // Not all values are used for each reaction type; see Table F.2 of MCNP User
  // Manual.
  std::vector<int> JXS;
  Datatable ESZ;   // energy table
  Datatable LONE;  // first word of table
  Datatable ITIE;  // inelastic energy table
  Datatable ESZG;  // energy table
  Datatable NU;    // fission nu data
  Datatable ITIX;  // inelastic cross sections
  Datatable JINC;  // incoherent form factors
  Datatable MTR;   // MT array
  Datatable ITXE;  // inelastic energy/angle distributions
  Datatable JCOH;  // coherent form factors
  Datatable LQR;   // Q-value array
  Datatable ITCE;  // elastic energy table
  Datatable JFLO;  // fluorescence data
  Datatable TYR;   // reaction type array
  Datatable ITCX;  // elastic cross sections
  Datatable LHNM;  // heating numbers
  Datatable LSIG;  // cross-section locators
  Datatable ITCA;  // elastic angular distributions
  Datatable LNEPS; // the number of electrons per shell
  Datatable SIG;   // cross sections
  Datatable SIGD;  // cross sections
  Datatable LBEPS; // binding energy per shell
  Datatable LAND;  // table of angular distribution locators
  Datatable LPIPS; // probability of interaction per shell
  Datatable AND;   // angular distributions
  Datatable LSWD;  // array of offsets to shellwise data
  Datatable LDLW;  // table of energy distribution locators
  Datatable SWD;   // shellwise data in PDF and CDF form
  Datatable DLW;   // energy distributions
  Datatable GPD;   // photon production data
  Datatable MTRP;  // photon production MT array
  Datatable LSIGP; // table of photon production cross-section locators
  Datatable SIGP;  // photon production cross sections
  Datatable LANDP; // table of photon production angular distribution
                   // locators
  Datatable ANDP;  // photon production angular distributions
  Datatable LDLWP; // table of photon production energy distribution
                   // locators
  Datatable DLWP;  // photon production energy distributions
  Datatable YP;    // table of yield multipliers
  Datatable FIS;   // total fission cross section
  Datatable END;   // last word of this table
  Datatable LUNR;  // probability tables
  Datatable DNU;   // delayed nubar data
  Datatable BDD;   // basic delayed data (lambdas, probabilities)
  Datatable DNEDL; // table of energy distribution locators
  Datatable DNED;  // energy distributions
};
} // end namespace rtt_nuclear_data

#endif // __nuclear_data_Nuclear_Data_hh__

//---------------------------------------------------------------------------//
// end of nuclear_data/Nuclear_Data.hh
//---------------------------------------------------------------------------//
