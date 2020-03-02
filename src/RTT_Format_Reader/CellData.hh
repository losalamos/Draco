//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   RTT_Format_Reader/CellData.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/CellData class.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_CellData_hh
#define rtt_RTT_Format_Reader_CellData_hh

#include "Cells.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the
 *        celldata block of the mesh file.
 */
class CellData {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;
  typedef std::vector<double> vector_dbl;
  typedef std::vector<std::vector<double>> vector_vector_dbl;

  const Dims &dims;
  vector_vector_dbl data;

public:
  CellData(const Dims &dims_)
      : dims(dims_),
        data(dims.get_ncells(), vector_dbl(dims.get_ncell_data())) {}
  ~CellData() {}

  void readCellData(ifstream &meshfile);

private:
  void readKeyword(ifstream &meshfile);
  void readData(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
   * \brief Returns all of the data field values for each of the cells.
   * \return The data field values for each of the cells.
   */
  vector_vector_dbl get_data() const { return data; }

  /*!
   * \brief Returns all of the data field values for the specified cell.
   * \param cell_numb Cell number.
   * \return The cell data field values.
   */
  vector_dbl get_data(size_t cell_numb) const { return data[cell_numb]; }

  /*!
   * \brief Returns the specified data field value for the specified cel.
   * \param cell_numb Cell number.
   * \param data_index Data field.
   * \return The cell data field value.
   */
  double get_data(size_t cell_numb, size_t data_index) const {
    return data[cell_numb][data_index];
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_CellData_hh

//----------------------------------------------------------------------------//
// end of RTT_Format_Reader/CellData.hh
//----------------------------------------------------------------------------//
