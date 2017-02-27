//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/Dims.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Dims class.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_Dims_hh__
#define __RTT_Format_Reader_Dims_hh__

#include "ds++/Assert.hh"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Controls parsing, storing, and accessing the data contained in the
 *        dimensions block of the mesh file.
 */
class Dims {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;
  typedef std::vector<int> vector_int;

  string coor_units;
  string prob_time_units;

  int ncell_defs;
  int nnodes_max;
  int nsides_max;
  int nnodes_side_max;

  int ndim;
  int ndim_topo;

  int nnodes;
  int nnode_flag_types;
  vector_int nnode_flags;
  int nnode_data;

  int nsides;
  int nside_types;
  vector_int side_types;
  int nside_flag_types;
  vector_int nside_flags;
  int nside_data;

  int ncells;
  int ncell_types;
  vector_int cell_types;
  int ncell_flag_types;
  vector_int ncell_flags;
  int ncell_data;

public:
  Dims(void)
      : coor_units(std::string()), prob_time_units(std::string()),
        ncell_defs(0), nnodes_max(0), nsides_max(0), nnodes_side_max(0),
        ndim(0), ndim_topo(0), nnodes(0), nnode_flag_types(0),
        nnode_flags(std::vector<int>()), nnode_data(0), nsides(0),
        nside_types(0), side_types(std::vector<int>()), nside_flag_types(0),
        nside_flags(std::vector<int>()), nside_data(0), ncells(0),
        ncell_types(0), cell_types(std::vector<int>()), ncell_flag_types(0),
        ncell_flags(std::vector<int>()), ncell_data(0) { /* empty */
  }
  ~Dims(void) { /*empty*/
  }

  void readDims(ifstream &meshfile);

private:
  void readKeyword(ifstream &meshfile);
  void readUnits(ifstream &meshfile);
  void readCellDefs(ifstream &meshfile);
  void readDimensions(ifstream &meshfile);
  void readNodes(ifstream &meshfile);
  void readSides(ifstream &meshfile);
  void readCells(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  // units and cell definition data access
  /*!
 * \brief Returns the problem coordinate units (e.g, cm).
 * \return Coordinate units.
 */
  string get_coor_units() const { return coor_units; }
  /*!
 * \brief Returns the problem time units (e.g, shakes).
 * \return Time units.
 */
  string get_prob_time_units() const { return prob_time_units; }
  /*!
 * \brief Returns the number of unique cell type definitions.
 * \return The number of cell definitions.
 */
  int get_ncell_defs() const { return ncell_defs; }
  /*!
 * \brief Returns the maximum number of nodes per cell type.
 * \return The maximum number of nodes per cell type.
 */
  int get_nnodes_max() const { return nnodes_max; }
  /*!
 * \brief Returns the maximum number of sides per cell type.
 * \return The maximum number of sides per cell type.
 */
  int get_nsides_max() const { return nsides_max; }
  /*!
 * \brief Returns the maximum number of nodes per cell side.
 * \return The maximum number of nodes per cell side.
 */
  int get_nnodes_side_max() const { return nnodes_side_max; }

  // dimensions node data access
  /*!
 * \brief Returns the number of spatial dimensions.
 * \return The number of spatial dimensions.
 */
  int get_ndim() const { return ndim; }
  /*!
 * \brief Returns the number of topological dimensions.
 * \return The number of topological dimensions.
 */
  int get_ndim_topo() const { return ndim_topo; }
  /*!
 * \brief Returns the number of nodes.
 * \return The number of nodes.
 */
  int get_nnodes() const { return nnodes; }
  /*!
 * \brief Returns the number of node flag types.
 * \return The number of node flag types.
 */
  int get_nnode_flag_types() const { return nnode_flag_types; }
  /*!
 * \brief Returns the number of node flags for the specified node flag type.
 * \param i Node flag type number.
 * \return The number of node flags.
 */
  int get_nnode_flags(int i) const { return nnode_flags[i]; }
  /*!
 * \brief Returns the number of node data fields.
 * \return The number of node data fields.
 */
  int get_nnode_data() const { return nnode_data; }

  // dimensions side data access
  /*!
 * \brief Returns the number of sides read.
 * \return The number of sides.
 */
  int get_nsides() const { return nsides; }
  /*!
 * \brief Returns the number of side types that are present in the "sides"
 *        block.
 * \return The number of side types.
 */
  int get_nside_types() const { return nside_types; }
  /*!
 * \brief Returns the side type index for the specified side type.
 * \param i Side type number.
 * \return The side type index.
 */
  int get_side_types(int i) const { return side_types[i]; }
  /*!
 * \brief Returns the number of side flag types.
 * \return The number of side flag types.
 */
  int get_nside_flag_types() const { return nside_flag_types; }
  /*!
 * \brief Returns the number of side flags for the specified side flag type.
 * \param i Side flag type number.
 * \return The number of side flags.
 */
  int get_nside_flags(int i) const { return nside_flags[i]; }
  /*!
 * \brief Returns the number of side data fields.
 * \return The number of side data fields.
 */
  int get_nside_data() const { return nside_data; }

  // cell data access
  /*!
 * \brief Returns the number of cells.
 * \return The number of cells.
 */
  int get_ncells() const { return ncells; }
  /*!
 * \brief Returns the number of cell types that are present in the "cells"
 *        block.
 * \return The number of cell types.
 */
  int get_ncell_types() const { return ncell_types; }
  /*!
 * \brief Returns the cell type index for the specified cell type.
 * \param i Cell type number.
 * \return The cell type index.
 */
  int get_cell_types(int i) const { return cell_types[i]; }
  /*!
 * \brief Returns the number of cell flag types.
 * \return The number of cell flag types.
 */
  int get_ncell_flag_types() const { return ncell_flag_types; }
  /*!
 * \brief Returns the number of cell flags for the specified cell flag type.
 * \param i Cell flag type number.
 * \return The number of cell flags.
 */
  int get_ncell_flags(int i) const { return ncell_flags[i]; }
  /*!
 * \brief Returns the number of cell data fields.
 * \return The number of cell data fields.
 */
  int get_ncell_data() const { return ncell_data; }

  /*!
 * \brief Validates the specified side type.
 * \param sidetype Side type number.
 * \return The existance of the side type.
 */
  bool allowed_side_type(int sidetype) const {
    return side_types.end() !=
           std::find(side_types.begin(), side_types.end(), sidetype);
  }
  /*!
 * \brief Validates the specified cell type.
 * \param celltype Cell type number.
 * \return The existance of the cell type.
 */
  bool allowed_cell_type(int celltype) const {
    return cell_types.end() !=
           std::find(cell_types.begin(), cell_types.end(), celltype);
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // __RTT_Format_Reader_Dims_hh__

//---------------------------------------------------------------------------//
// end of RTT_Format_Reader/Dims.hh
//---------------------------------------------------------------------------//
