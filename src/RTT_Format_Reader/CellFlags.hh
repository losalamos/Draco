//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/CellFlags.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/CellFlags class.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_CellFlags_hh__
#define __RTT_Format_Reader_CellFlags_hh__

#include "Dims.hh"
#include "Flags.hh"
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the
 *        cell flags block of the mesh file.
 */
class CellFlags {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;

  const Dims &dims;
  std::vector<std::shared_ptr<Flags>> flagTypes;

public:
  CellFlags(const Dims &dims_)
      : dims(dims_), flagTypes(dims.get_ncell_flag_types()) {}
  ~CellFlags() {}

  void readCellFlags(ifstream &meshfile);

private:
  void readKeyword(ifstream &meshfile);
  void readFlagTypes(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
 * \brief Validates the specified cell flag type and number.
 * \param flagtype Cell flag type number.
 * \param flag Flag number.
 * \return The existance of the cell flag type and number.
 */
  bool allowed_flag(int flagtype, int flag) const {
    Insist(flagtype <= dims.get_ncell_flag_types() - 1,
           "Invalid cell flag type number!");
    return flagTypes[flagtype]->allowed_flag(flag);
  }
  /*!
 * \brief Returns the name of specified cell flag type read from the mesh file
 *        cell_flags data.
 * \param flagtype Cell flag type number.
 * \return The cell flag type name.
 */
  string get_flag_type(int flagtype) const {
    Insist(flagtype <= dims.get_ncell_flag_types() - 1,
           "Invalid cell flag type number!");
    return flagTypes[flagtype]->getFlagType();
  }

  int get_flag_type_index(string &desired_flag_type) const;
  /*!
 * \brief Returns the cell flag number associated with the specified cell flag
 *        type and cell flag index.
 * \param flagtype Cell flag type number.
 * \param flag_index Cell flag index.
 * \return The cell flag number.
 */
  int get_flag_number(int flagtype, int flag_index) const {
    Insist(flagtype <= dims.get_ncell_flag_types() - 1,
           "Invalid cell flag type number!");
    Insist(flag_index <= flagTypes[flagtype]->getFlagSize() - 1,
           "Invalid cell flag number index number!");
    return flagTypes[flagtype]->getFlagNumber(flag_index);
  }
  /*!
 * \brief Returns the number of cell flags for the specified cell flag type.
 * \param flagtype Cell flag type number.
 * \return The number of cell flags.
 */
  int get_flag_size(int flagtype) const {
    Insist(flagtype <= dims.get_ncell_flag_types() - 1,
           "Invalid cell flag type number!");
    return flagTypes[flagtype]->getFlagSize();
  }
  /*!
 * \brief Returns the cell flag name associated with the specified cell flag
 *        type and cell flag index.
 * \param flagtype Cell flag type number.
 * \param flag_index Cell flag index.
 * \return The cell flag name.
 */
  string get_flag_name(int flagtype, int flag_index) const {
    Insist(flagtype <= dims.get_ncell_flag_types() - 1,
           "Invalid cell flag type number!");
    Insist(flag_index <= flagTypes[flagtype]->getFlagSize() - 1,
           "Invalid cell flag name index number!");
    return flagTypes[flagtype]->getFlagName(flag_index);
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // __RTT_Format_Reader_CellFlags_hh__

//---------------------------------------------------------------------------//
// end of RTT_Format_Reader/CellFlags.hh
//---------------------------------------------------------------------------//
