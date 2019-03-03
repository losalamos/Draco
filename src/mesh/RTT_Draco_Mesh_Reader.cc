//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/RTT_Draco_Mesh_Reader.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Friday, Jul 13, 2018, 08:38 am
 * \brief  RTT_Draco_Mesh_Reader header file.
 * \note   Copyright (C) 2018-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "RTT_Draco_Mesh_Reader.hh"

namespace rtt_mesh {

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief RTT_Draco_Mesh_Reader constructor.
 *
 * \param[in] filename_ name of file to be parsed
 */
RTT_Draco_Mesh_Reader::RTT_Draco_Mesh_Reader(const std::string filename_)
    : filename(filename_) {
  // check for valid file name
  Insist(filename_.size() > 0, "No file name supplied.");

  // \todo: remove read_mesh function and read in constructor(?)
  Require(rtt_reader == nullptr);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Read the mesh by constructing an RTT_Format_Reader
 */
void RTT_Draco_Mesh_Reader::read_mesh() {

  // insist the file can be opened before trying read_mesh() wrapper
  std::ifstream rttfile(filename.c_str());
  Insist(rttfile.is_open(), "Failed to find or open specified RTT mesh file.");
  rttfile.close();

  rtt_reader.reset(new rtt_RTT_Format_Reader::RTT_Format_Reader(filename));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of nodes for a cell.
 *
 * \param[in] cell index of cell (0-based)
 *
 * \return number of nodes for cell
 */
unsigned RTT_Draco_Mesh_Reader::get_celltype(size_t cell) const {

  // first obtain a cell definition index
  size_t cell_def = rtt_reader->get_cells_type(cell);

  // for Draco_Mesh, cell_type is number of nodes
  Check(rtt_reader->get_cell_defs_nnodes(cell_def) < UINT_MAX);
  unsigned cell_type =
      static_cast<unsigned>(rtt_reader->get_cell_defs_nnodes(cell_def));

  return cell_type;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the vector of node indices for a side.
 *
 * \param[in] side index of side (0-based)
 *
 * \return number of nodes for side
 */
size_t RTT_Draco_Mesh_Reader::get_sidetype(size_t side) const {

  // first obtain a side definition index
  size_t side_def = rtt_reader->get_sides_type(side);

  // acquire the number of nodes associated with this side def
  size_t side_type = rtt_reader->get_cell_defs_nnodes(side_def);

  return side_type;
}

} // end namespace rtt_mesh

//---------------------------------------------------------------------------//
// end of mesh/RTT_Draco_Mesh_Reader.cc
//---------------------------------------------------------------------------//
