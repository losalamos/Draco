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

  const char *test = filepath.c_str();
  ifstream ACEfile(test, std::ios::in);
  if (!ACEfile) {
    std::ostringstream buffer;
    buffer << "File " << filepath << " could not be opened\n";
    throw std::invalid_argument(buffer.str());
  }

  string line;
  std::getline(ACEfile, line);
  zaid = line.substr(0, 10);
  atomic_weight = stod(line.substr(10, 12));
  temperature = stod(line.substr(22, 12));
  date = line.substr(35, 10); // Starting index may be wrong, manual is
                              // ambiguous.

  std::getline(ACEfile, line);
  comment = line.substr(0, 70);
  mat_identifier = line.substr(70, 10);

  // Reaction type given by char in ZAID
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

  atomic_number = stoi(zaid.substr(0, 3));
  mass_number = stoi(zaid.substr(3, 3)); // 000 for P, M, G, E
  evaluation_identifier = stoi(zaid.substr(7, 2));

  for (int n = 0; n < 4; n++) {
    std::getline(ACEfile, line);
    for (int m = 0; m < 4; m++) {
      IZ.push_back(stoi(line.substr(0 + 18 * m, 7)));
      AW.push_back(stod(line.substr(7 + 18 * m, 11)));
    }
  }

  string NXS_s, JXS_s, data_s;
  for (int n = 0; n < 2; n++) {
    std::getline(ACEfile, line);
    NXS_s.append(line);
  }
  for (int n = 0; n < 4; n++) {
    std::getline(ACEfile, line);
    JXS_s.append(line);
  }
  for (int n = 0; n < 16; n++) {
    NXS.push_back(stoi(NXS_s.substr(0 + 9 * n, 9)));
  }
  for (int n = 0; n < 32; n++) {
    JXS.push_back(stoi(JXS_s.substr(0 + 9 * n, 9)));
  }

  while (true) {
    std::getline(ACEfile, line);
    data_s.append(line);
    if (ACEfile.eof())
      break;
  }
  ACEfile.close();

  // Interpret file depending on reaction type
  if (reaction == Reaction::CONTINUOUS_ENERGY_NEUTRON ||
      reaction == Reaction::DISCRETE_REACTION_NEUTRON) {
    std::ostringstream buffer;
    buffer << "Reaction " << Reaction_Name[reaction] << " not supported\n";
    throw std::invalid_argument(buffer.str());
  } else if (reaction == Reaction::DOSIMETRY) {
    std::ostringstream buffer;
    buffer << "Reaction " << Reaction_Name[reaction] << " not supported\n";
    throw std::invalid_argument(buffer.str());
  } else if (reaction == Reaction::THERMAL_S) {
    std::ostringstream buffer;
    buffer << "Reaction " << Reaction_Name[reaction] << " not supported\n";
    throw std::invalid_argument(buffer.str());
  } else if (reaction == Reaction::CONTINUOUS_ENERGY_PHOTOATOMIC) {
    data_length = NXS[0];
    Z = NXS[1];
    NES = NXS[2];
    NFLO = NXS[3];
    NSH = NXS[4];
    NT = NXS[14];
    photon_production = NXS[15];

    ESZG.set_indices(JXS[0], JXS[1] - JXS[0]);
    JINC.set_indices(JXS[1], JXS[2] - JXS[1]);
    JCOH.set_indices(JXS[2], JXS[3] - JXS[2]);
    JFLO.set_indices(JXS[3], JXS[4] - JXS[3]);
    LHNM.set_indices(JXS[4], JXS[5] - JXS[4]);
    LNEPS.set_indices(JXS[5], JXS[6] - JXS[5]);
    LBEPS.set_indices(JXS[6], JXS[7] - JXS[6]);
    LPIPS.set_indices(JXS[7], JXS[8] - JXS[7]);
    LSWD.set_indices(JXS[8], JXS[9] - JXS[8]);
    SWD.set_indices(JXS[9], data_length - JXS[9]);
  } else if (reaction == Reaction::NEUTRON_MULTIGROUP) {
    std::ostringstream buffer;
    buffer << "Reaction " << Reaction_Name[reaction] << " not supported\n";
    throw std::invalid_argument(buffer.str());
  } else if (reaction == Reaction::PHOTOATOMIC_MULTIGROUP) {
    std::ostringstream buffer;
    buffer << "Reaction " << Reaction_Name[reaction] << " not supported\n";
    throw std::invalid_argument(buffer.str());
  } else if (reaction == Reaction::CONTINUOUS_ENERGY_ELECTRON) {
    std::ostringstream buffer;
    buffer << "Reaction " << Reaction_Name[reaction] << " not supported\n";
    throw std::invalid_argument(buffer.str());
  } else if (reaction == Reaction::CONTINUOUS_ENERGY_PHOTONUCLEAR) {
    std::ostringstream buffer;
    buffer << "Reaction " << Reaction_Name[reaction] << " not supported\n";
    throw std::invalid_argument(buffer.str());
  } else {
    std::ostringstream buffer;
    buffer << "Reaction " << Reaction_Name[reaction] << " not supported\n";
    throw std::invalid_argument(buffer.str());
  }

  // Give every datatable the chance to be read from file.
  load_datatable(ESZ, data_s);
  load_datatable(LONE, data_s);
  load_datatable(ITIE, data_s);
  load_datatable(ESZG, data_s);
  load_datatable(NU, data_s);
  load_datatable(ITIX, data_s);
  load_datatable(JINC, data_s);
  load_datatable(MTR, data_s);
  load_datatable(ITXE, data_s);
  load_datatable(JCOH, data_s);
  load_datatable(LQR, data_s);
  load_datatable(ITCE, data_s);
  load_datatable(JFLO, data_s);
  load_datatable(TYR, data_s);
  load_datatable(ITCX, data_s);
  load_datatable(LHNM, data_s);
  load_datatable(LSIG, data_s);
  load_datatable(ITCA, data_s);
  load_datatable(LNEPS, data_s);
  load_datatable(SIG, data_s);
  load_datatable(SIGD, data_s);
  load_datatable(LBEPS, data_s);
  load_datatable(LAND, data_s);
  load_datatable(LPIPS, data_s);
  load_datatable(AND, data_s);
  load_datatable(LSWD, data_s);
  load_datatable(LDLW, data_s);
  load_datatable(SWD, data_s);
  load_datatable(DLW, data_s);
  load_datatable(GPD, data_s);
  load_datatable(MTRP, data_s);
  load_datatable(LSIGP, data_s);
  load_datatable(SIGP, data_s);
  load_datatable(LANDP, data_s);
  load_datatable(ANDP, data_s);
  load_datatable(LDLWP, data_s);
  load_datatable(DLWP, data_s);
  load_datatable(YP, data_s);
  load_datatable(FIS, data_s);
  load_datatable(END, data_s);
  load_datatable(LUNR, data_s);
  load_datatable(DNU, data_s);
  load_datatable(BDD, data_s);
  load_datatable(DNEDL, data_s);
  load_datatable(DNED, data_s);
}

void Nuclear_Data::load_datatable(Datatable &dataTable,
                                  const std::string data_s) {
  // Only load data requested for this reaction type
  if (dataTable.start_index >= 0) {
    string subdata_s =
        data_s.substr(20 * dataTable.start_index, 20 * dataTable.length);
    int nentries = subdata_s.length() / 20;
    for (int n = 0; n < nentries; n++) {
      dataTable.data.push_back(stod(subdata_s.substr(20 * n, 20)));
    }
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
  std::cout << "Reaction:              " << Reaction_Name[reaction]
            << std::endl;
  std::cout << "IZ, AW pairs:          " << std::endl;
  for (int n = 0; n < 16; n++) {
    std::cout << "  " << IZ[n] << " " << AW[n] << std::endl;
  }
  std::cout << "NXS" << std::endl;
  std::cout << "data_length: " << data_length << std::endl;
}
} // end namespace rtt_nuclear_data
