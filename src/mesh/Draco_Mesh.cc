//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Draco_Mesh.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jun 07, 2018, 15:38 pm
 * \brief  Draco_Mesh class implementation file.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Draco_Mesh.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <iostream>
#include <numeric>

namespace rtt_mesh {

// helper function for safe conversion of types during initialization list
// processing.
unsigned safe_convert_from_size_t(size_t const in_) {
  Check(in_ < UINT_MAX);
  return static_cast<unsigned>(in_);
}

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Draco_Mesh constructor.
 *
 * \param[in] dimension_ dimension of mesh
 * \param[in] geometry_ enumerator of possible coordinate system geometries
 * \param[in] cell_type_ number of vertices for each cell
 * \param[in] cell_to_node_linkage_ serialized map of cell indices to node
 *               indices.
 * \param[in] side_set_flag_ map of side indices (per cell) to side flag (global
 *               index for a side).
 * \param[in] side_node_count_ number of nodes per each cell on a side of
 *               the mesh.
 * \param[in] side_to_node_linkage_ serialized map of side indices (per side
 *               cell) to node indices.
 * \param[in] coordinates_ serialized map of node index to coordinate values.
 * \param[in] global_node_number_ map of local to global node index (vector
 *               subscript is local node index and value is global node index;
 *               for one process, this is the identity map).
 * \param[in] ghost_cell_type_ number of vertices corresponding to each ghost
 *               cell (1 in 1D, 2 in 2D, arbitrary in 3D).
 * \param[in] ghost_cell_to_node_linkage_ serialized map of index into vector of
 *               ghost cells to local index of ghost nodes.
 * \param[in] ghost_cell_number_ cell index local to other processor.
 * \param[in] ghost_cell_rank_ rank of each ghost cell.
 */
Draco_Mesh::Draco_Mesh(unsigned dimension_, Geometry geometry_,
                       const std::vector<unsigned> &cell_type_,
                       const std::vector<unsigned> &cell_to_node_linkage_,
                       const std::vector<unsigned> &side_set_flag_,
                       const std::vector<unsigned> &side_node_count_,
                       const std::vector<unsigned> &side_to_node_linkage_,
                       const std::vector<double> &coordinates_,
                       const std::vector<unsigned> &global_node_number_,
                       const std::vector<unsigned> &ghost_cell_type_,
                       const std::vector<unsigned> &ghost_cell_to_node_linkage_,
                       const std::vector<unsigned> &ghost_cell_number_,
                       const std::vector<unsigned> &ghost_cell_rank_)
    : dimension(dimension_), geometry(geometry_),
      num_cells(safe_convert_from_size_t(cell_type_.size())),
      num_nodes(safe_convert_from_size_t(global_node_number_.size())),
      side_set_flag(side_set_flag_), ghost_cell_number(ghost_cell_number_),
      ghost_cell_rank(ghost_cell_rank_),
      node_coord_vec(compute_node_coord_vec(coordinates_)),
      m_cell_type(cell_type_), m_cell_to_node_linkage(cell_to_node_linkage_),
      m_side_node_count(side_node_count_),
      m_side_to_node_linkage(side_to_node_linkage_) {

  // Require(dimension_ <= 3);
  // \todo: generalize mesh generation to 1D,3D (and uncomment requirment above)
  Insist(dimension_ == 2, "dimension_ != 2");

  // require some constraints on vector sizes
  Require(cell_to_node_linkage_.size() ==
          std::accumulate(cell_type_.begin(), cell_type_.end(), 0u));
  Require(
      side_to_node_linkage_.size() ==
      std::accumulate(side_node_count_.begin(), side_node_count_.end(), 0u));
  Require(coordinates_.size() == dimension_ * global_node_number_.size());

  // check ghost data (should be true even when none are supplied)
  Require(ghost_cell_type_.size() == ghost_cell_number_.size());
  Require(ghost_cell_rank_.size() == ghost_cell_number_.size());
  Require(
      ghost_cell_to_node_linkage_.size() ==
      std::accumulate(ghost_cell_type_.begin(), ghost_cell_type_.end(), 0u));

  // build the layout (or linkage) of the mesh
  compute_cell_to_cell_linkage(cell_type_, cell_to_node_linkage_,
                               side_node_count_, side_to_node_linkage_,
                               ghost_cell_type_, ghost_cell_to_node_linkage_);
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 *
 * \param[in] coordinates serialized map of node index to coordinate values
 *               (passed from constructor).
 * \return a vector of vectors of size=dimension of coordinates.
 */
std::vector<std::vector<double>> Draco_Mesh::compute_node_coord_vec(
    const std::vector<double> &coordinates) const {

  Require(coordinates.size() == dimension * num_nodes);

  // resize this class's coordinate data member
  std::vector<std::vector<double>> ret_node_coord_vec(
      num_nodes, std::vector<double>(dimension));

  // de-serialize the vector of node coordinates
  std::vector<double>::const_iterator ncv_first = coordinates.begin();
  for (unsigned node = 0; node < num_nodes; ++node) {

    // use the cell_type to create a vector of node indices for this cell
    std::vector<double> coord_vec(ncv_first, ncv_first + dimension);

    // resize each entry to the number of dimensions
    ret_node_coord_vec[node] = coord_vec;

    // increment pointer
    ncv_first += dimension;
  }

  Ensure(ncv_first == coordinates.end());

  return ret_node_coord_vec;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 *
 * \param[in] cell_type number of vertices per cell.
 * \param[in] cell_to_node_linkage serial map of cell index to node indices.
 * \param[in] side_node_count number of verices per side.
 * \param[in] side_to_node_linkage serial map of side index to node indices.
 */
void Draco_Mesh::compute_cell_to_cell_linkage(
    const std::vector<unsigned> &cell_type,
    const std::vector<unsigned> &cell_to_node_linkage,
    const std::vector<unsigned> &side_node_count,
    const std::vector<unsigned> &side_to_node_linkage,
    const std::vector<unsigned> &ghost_cell_type,
    const std::vector<unsigned> &ghost_cell_to_node_linkage) {

  Require(cell_type.size() == num_cells);
  Require(cell_to_node_linkage.size() ==
          std::accumulate(cell_type.begin(), cell_type.end(), 0u));
  Require(side_to_node_linkage.size() ==
          std::accumulate(side_node_count.begin(), side_node_count.end(), 0u));

  // STEP 1: create a node-to-cell map (inverse of step 1)

  std::map<unsigned, std::vector<unsigned>> node_to_cell_map =
      compute_node_indx_map(cell_type, cell_to_node_linkage);

  // STEP 2: create a node-to-side map

  std::map<unsigned, std::vector<unsigned>> node_to_side_map =
      compute_node_indx_map(side_node_count, side_to_node_linkage);

  Remember(const size_t num_sides =
               safe_convert_from_size_t(side_node_count.size()));
  Check(dimension == 2 ? side_to_node_linkage.size() == 2 * num_sides : true);
  Check(node_to_side_map.size() == num_sides);

  // STEP 3: create a node-to-ghost-cell map

  std::map<unsigned, std::vector<unsigned>> node_to_ghost_cell_map =
      compute_node_indx_map(ghost_cell_type, ghost_cell_to_node_linkage);

  // STEP 4: identify faces and create cell to cell map

  // \todo: global face index?
  // \todo: extend to 1D, 3D
  // \todo: add all cells as keys to each layout, even without values

  // identify faces per cell and create cell-to-cell linkage
  // in 2D, faces will always have 2 nodes
  unsigned node_offset = 0;
  for (unsigned cell = 0; cell < num_cells; ++cell) {

    // create a vector of all possible node pairs
    unsigned nper_cell = cell_type[cell];
    unsigned num_pairs = nper_cell * (nper_cell - 1) / 2;
    std::vector<std::vector<unsigned>> vec_node_vec(num_pairs,
                                                    std::vector<unsigned>(2));
    // \todo: change type of vec_node_vec?
    unsigned k = 0;
    for (unsigned i = 0; i < nper_cell; ++i) {
      for (unsigned j = i + 1; j < nper_cell; ++j) {

        // set kth pair entry to the size 2 vector of nodes
        vec_node_vec[k] = {cell_to_node_linkage[node_offset + i],
                           cell_to_node_linkage[node_offset + j]};

        // increment k
        k++;
      }
    }

    Check(k == num_pairs);

    // check if each pair constitutes a face
    for (unsigned l = 0; l < num_pairs; ++l) {

      // get adjacent cells from node-to-cell map
      // \todo: add DbC to ensure these are sorted
      const std::vector<unsigned> &vert0_cells =
          node_to_cell_map[vec_node_vec[l][0]];
      const std::vector<unsigned> &vert1_cells =
          node_to_cell_map[vec_node_vec[l][1]];

      // find common cells (set_intersection is low-complexity)
      // \todo: reserve size for cells_in_common
      std::vector<unsigned> cells_in_common;
      std::set_intersection(vert0_cells.begin(), vert0_cells.end(),
                            vert1_cells.begin(), vert1_cells.end(),
                            std::back_inserter(cells_in_common));

      // these nodes should have at least cell index "cell" in common
      Check(cells_in_common.size() >= 1);

      // populate cell-to-cell linkage
      if (cells_in_common.size() > 1) {
        for (auto oth_cell : cells_in_common) {
          if (oth_cell != cell)
            cell_to_cell_linkage[cell].push_back(
                std::make_pair(oth_cell, vec_node_vec[l]));
        }
      }

      // check if this vertex pair has a side flag
      // \todo: add DbC to ensure these are sorted
      const std::vector<unsigned> &vert0_sides =
          node_to_side_map[vec_node_vec[l][0]];
      const std::vector<unsigned> &vert1_sides =
          node_to_side_map[vec_node_vec[l][1]];

      // find common sides
      std::vector<unsigned> sides_in_common;
      std::set_intersection(vert0_sides.begin(), vert0_sides.end(),
                            vert1_sides.begin(), vert1_sides.end(),
                            std::back_inserter(sides_in_common));

      Check(sides_in_common.size() <= 1);
      if (sides_in_common.size() > 0) {
        // populate cell-to-side linkage
        cell_to_side_linkage[cell].push_back(
            std::make_pair(sides_in_common[0], vec_node_vec[l]));
      }

      // check if this vertex pair has a ghost cell
      // \todo: add DbC to ensure these are sorted
      const std::vector<unsigned> &vert0_ghosts =
          node_to_ghost_cell_map[vec_node_vec[l][0]];
      const std::vector<unsigned> &vert1_ghosts =
          node_to_ghost_cell_map[vec_node_vec[l][1]];

      // find common ghost cells
      std::vector<unsigned> ghost_cells_in_common;
      std::set_intersection(vert0_ghosts.begin(), vert0_ghosts.end(),
                            vert1_ghosts.begin(), vert1_ghosts.end(),
                            std::back_inserter(ghost_cells_in_common));

      Check(ghost_cells_in_common.size() <= 1);
      if (ghost_cells_in_common.size() > 0) {
        // populated cell-to-ghost-cell linkage
        cell_to_ghost_cell_linkage[cell].push_back(
            std::make_pair(ghost_cells_in_common[0], vec_node_vec[l]));
      }
    }

    // increment serialized cell_to_node linkage index offset
    node_offset += nper_cell;
  }

  Check(node_offset == cell_to_node_linkage.size());

  // STEP 5: instantiate the full layout
  // \todo: finish Draco_Layout class

  Ensure(cell_to_cell_linkage.size() <= num_cells);
  Ensure(cell_to_side_linkage.size() <= num_cells);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build an intermediate node-index map to support layout generation.
 *
 * \param[in] indx_type vector of number of nodes, subscipted by index.
 * \param[in] indx_to_node_linkage serial map of index to node indices.
 * \return a map of node index to vector of indexes adjacent to the node.
 */
std::map<unsigned, std::vector<unsigned>> Draco_Mesh::compute_node_indx_map(
    const std::vector<unsigned> &indx_type,
    const std::vector<unsigned> &indx_to_node_linkage) const {

  // map to return
  std::map<unsigned, std::vector<unsigned>> node_to_indx_map;

  // push node index to vector for each node
  unsigned node_offset = 0;
  const size_t num_indxs = indx_type.size();
  for (unsigned indx = 0; indx < num_indxs; ++indx) {

    // push the indx onto indx vectors for each node
    for (unsigned i = 0; i < indx_type[indx]; ++i) {

      Check(indx_to_node_linkage[node_offset + i] < num_nodes);

      node_to_indx_map[indx_to_node_linkage[node_offset + i]].push_back(indx);
    }

    // increment offset
    node_offset += indx_type[indx];
  }

  Ensure(node_offset == indx_to_node_linkage.size());

  return node_to_indx_map;
}

} // end namespace rtt_mesh

//---------------------------------------------------------------------------//
// end of mesh/Draco_Mesh.cc
//---------------------------------------------------------------------------//
