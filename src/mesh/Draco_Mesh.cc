//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   mesh/Draco_Mesh.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jun 07, 2018, 15:38 pm
 * \brief  Draco_Mesh class implementation file.
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Draco_Mesh.hh"
#include "c4/C4_Functions.hh"
#include "c4/gatherv.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <numeric>

namespace rtt_mesh {

// helper function for safe conversion of types during initialization list processing.
unsigned safe_convert_from_size_t(size_t const in_) {
  Check(in_ < UINT_MAX);
  return static_cast<unsigned>(in_);
}

//------------------------------------------------------------------------------------------------//
// CONSTRUCTOR
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Draco_Mesh constructor.
 *
 * \param[in] dimension_ dimension of mesh
 * \param[in] geometry_ enumerator of possible coordinate system geometries
 * \param[in] num_faces_per_cell_ number of faces for each cell.
 * \param[in] cell_to_node_linkage_ serialized map of cell indices to node indices. nodes are listed
 *               per cell face. so there are duplicate node entries in 2D or 3D since adjacent cell
 *               faces will share one or more nodes. 2D node ordering will be assumed to be
 *               counterclockwise around the cell, in 3D the node ordering per face is assumed to be
 *               counterclockwise from inside the cell looking at the face.
 * \param[in] side_set_flag_ map of side indices (per cell) to side flag (global index for a side).
 * \param[in] side_node_count_ number of nodes per each cell on a side of the mesh.
 * \param[in] side_to_node_linkage_ serialized map of side indices (per side cell) to node indices.
 * \param[in] coordinates_ serialized map of node index to coordinate values.
 * \param[in] global_node_number_ map of local to global node index (vector subscript is local node
 *               index and value is global node index; for one process, this is the identity map).
 * \param[in] num_nodes_per_face_per_cell_ number of vertices per face per cell.
 * \param[in] ghost_cell_type_ number of vertices corresponding to each ghost cell (1 in 1D, 2 in
 *               2D, arbitrary in 3D).
 * \param[in] ghost_cell_to_node_linkage_ serialized map of index into vector of ghost cells to
 *               local index of ghost nodes.
 * \param[in] ghost_cell_number_ cell index local to other processor.
 * \param[in] ghost_cell_rank_ rank of each ghost cell.
 */
Draco_Mesh::Draco_Mesh(
    unsigned dimension_, Geometry geometry_, const std::vector<unsigned> &num_faces_per_cell_,
    const std::vector<unsigned> &cell_to_node_linkage_, const std::vector<unsigned> side_set_flag_,
    const std::vector<unsigned> &side_node_count_,
    const std::vector<unsigned> &side_to_node_linkage_, const std::vector<double> &coordinates_,
    const std::vector<unsigned> &global_node_number_,
    const std::vector<unsigned> &num_nodes_per_face_per_cell_,
    const std::vector<unsigned> &ghost_cell_type_,
    const std::vector<unsigned> &ghost_cell_to_node_linkage_,
    const std::vector<int> &ghost_cell_number_, const std::vector<int> &ghost_cell_rank_)
    : dimension(dimension_), geometry(geometry_),
      num_cells(safe_convert_from_size_t(num_faces_per_cell_.size())),
      num_nodes(safe_convert_from_size_t(global_node_number_.size())),
      side_set_flag(std::move(side_set_flag_)), ghost_cell_number(ghost_cell_number_),
      ghost_cell_rank(ghost_cell_rank_), node_coord_vec(compute_node_coord_vec(coordinates_)),
      m_num_faces_per_cell(num_faces_per_cell_),
      m_num_nodes_per_face_per_cell(num_nodes_per_face_per_cell_),
      m_cell_to_node_linkage(compute_cell_to_node_tensor(
          num_faces_per_cell_, num_nodes_per_face_per_cell_, cell_to_node_linkage_)),
      m_side_node_count(side_node_count_), m_side_to_node_linkage(side_to_node_linkage_) {

  Require(dimension_ <= 3);
  Require(side_to_node_linkage_.size() ==
          std::accumulate(side_node_count_.begin(), side_node_count_.end(), 0u));
  Require(coordinates_.size() == dimension_ * global_node_number_.size());

  // check ghost data (should be true even when none are supplied)
  Require(ghost_cell_type_.size() == ghost_cell_number_.size());
  Require(ghost_cell_rank_.size() == ghost_cell_number_.size());
  Require(ghost_cell_to_node_linkage_.size() ==
          std::accumulate(ghost_cell_type_.begin(), ghost_cell_type_.end(), 0u));

  // build the layout using face types (number of nodes per face per cell)
  compute_cell_to_cell_linkage(
      num_faces_per_cell_, cell_to_node_linkage_, num_nodes_per_face_per_cell_, side_node_count_,
      side_to_node_linkage_, ghost_cell_type_, ghost_cell_to_node_linkage_);

  // build cell-to-corner-cell layout
  compute_node_to_cell_linkage(num_faces_per_cell_, cell_to_node_linkage_,
                               num_nodes_per_face_per_cell_, ghost_cell_type_,
                               ghost_cell_to_node_linkage_, global_node_number_);
}

//------------------------------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Obtain a unique list of a cell's nodes.
 *
 * \param[in] cell the index of the cell
 *
 * \return a vector of node indices for the cell, without duplicates
 */
const std::vector<unsigned> Draco_Mesh::get_cell_nodes(const unsigned cell) const {
  Require(cell < num_cells);

  // initialize return vector
  std::vector<unsigned> ret_cell_nodes;

  for (unsigned face = 0; face < m_num_faces_per_cell[cell]; ++face) {

    for (auto node : m_cell_to_node_linkage[cell][face]) {

      // this preserves counter-clockwise ordering in 2D
      if (std::find(ret_cell_nodes.begin(), ret_cell_nodes.end(), node) == ret_cell_nodes.end())
        ret_cell_nodes.push_back(node);
    }
  }

  Ensure(!ret_cell_nodes.empty());
  return ret_cell_nodes;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Return a flattened version of the cell-node tensor
 *
 * \return a flattened cell-node linkage vector
 */
const std::vector<unsigned> Draco_Mesh::get_flat_cell_node_linkage() const {

  // initialize empty cell-node vector
  std::vector<unsigned> ret_flat_cell_node;

  // insert node vectors from each face of each cell
  const size_t cn_size = m_cell_to_node_linkage.size();
  for (size_t cell = 0; cell < cn_size; ++cell) {

    const size_t cnface_size = m_cell_to_node_linkage[cell].size();
    for (size_t face = 0; face < cnface_size; ++face) {

      const std::vector<unsigned> node_vec = m_cell_to_node_linkage[cell][face];
      ret_flat_cell_node.insert(ret_flat_cell_node.end(), node_vec.begin(), node_vec.end());
    }
  }

  // return the flattened cell-node tensor
  Ensure(ret_flat_cell_node.size() >= 2 * m_cell_to_node_linkage.size());
  Ensure(!ret_flat_cell_node.empty());
  return ret_flat_cell_node;
}

//------------------------------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 *
 * \param[in] coordinates serialized map of node index to coordinate values (passed from
 *               constructor).
 * \return a vector of vectors of size=dimension of coordinates.
 */
std::vector<std::vector<double>>
Draco_Mesh::compute_node_coord_vec(const std::vector<double> &coordinates) const {

  Require(coordinates.size() == dimension * static_cast<size_t>(num_nodes));

  // resize this class's coordinate data member
  std::vector<std::vector<double>> ret_node_coord_vec(num_nodes, std::vector<double>(dimension));

  // deserialize the vector of node coordinates
  auto ncv_first = coordinates.begin();
  for (unsigned node = 0; node < num_nodes; ++node) {

    // create a vector of node indices for this cell
    std::vector<double> coord_vec(ncv_first, ncv_first + dimension);

    // resize each entry to the number of dimensions
    ret_node_coord_vec[node] = coord_vec;

    // increment pointer
    ncv_first += dimension;
  }

  Ensure(ncv_first == coordinates.end());

  return ret_node_coord_vec;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 *
 * \param[in] num_faces_per_cell number of faces per cell
 * \param[in] num_nodes_per_face_per_cell number of nodes per face per cell.
 * \param[in] cell_to_node_linkage node indices per face per cell (serialized).
 *
 * \return a vector of vectors of size=dimension of coordinates.
 */
std::vector<std::vector<std::vector<unsigned>>>
Draco_Mesh::compute_cell_to_node_tensor(const std::vector<unsigned> &num_faces_per_cell,
                                        const std::vector<unsigned> &num_nodes_per_face_per_cell,
                                        const std::vector<unsigned> &cell_to_node_linkage) const {

  std::vector<std::vector<std::vector<unsigned>>> ret_cn_tensor(num_cells);

  // deserialize the cell-node linkage vector
  auto cfn_first = cell_to_node_linkage.begin();
  unsigned cf_indx = 0;
  for (unsigned cell = 0; cell < num_cells; ++cell) {

    // resize the vector of vectors of nodes to the number of faces
    ret_cn_tensor[cell].resize(num_faces_per_cell[cell]);

    for (unsigned face = 0; face < num_faces_per_cell[cell]; ++face) {

      // get the vector of nodes
      std::vector<unsigned> node_vec(cfn_first, cfn_first + num_nodes_per_face_per_cell[cf_indx]);

      // resize and set the vector of nodes for this face
      ret_cn_tensor[cell][face] = node_vec;

      // increment iterator and counter
      cfn_first += num_nodes_per_face_per_cell[cf_indx];
      cf_indx++;
    }
  }

  Ensure(cfn_first == cell_to_node_linkage.end());

  return ret_cn_tensor;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build the cell-face index map to the corresponding coordinates.
 *
 * \param[in] num_faces_per_cell number of faces per cell.
 * \param[in] cell_to_node_linkage serial map of cell to face to node indices.
 * \param[in] num_nodes_per_face_per_cell number of nodes per face per cell
 * \param[in] side_node_count number of vertices per side.
 * \param[in] side_to_node_linkage serial map of side index to node indices.
 * \param[in] ghost_cell_type  number of common vertices per ghost cell.
 * \param[in] ghost_cell_to_node_linkage vertices in common per ghost cell.
 */
void Draco_Mesh::compute_cell_to_cell_linkage(
    const std::vector<unsigned> &num_faces_per_cell,
    const std::vector<unsigned> &cell_to_node_linkage,
    const std::vector<unsigned> &num_nodes_per_face_per_cell,
    const std::vector<unsigned> &side_node_count, const std::vector<unsigned> &side_to_node_linkage,
    const std::vector<unsigned> &ghost_cell_type,
    const std::vector<unsigned> &ghost_cell_to_node_linkage) {

  Require(num_nodes_per_face_per_cell.size() > 0);
  Require(num_nodes_per_face_per_cell.size() ==
          std::accumulate(num_faces_per_cell.begin(), num_faces_per_cell.end(), 0u));

  // (1) create map of cell face to node set

  std::map<unsigned, std::set<unsigned>> cface_to_nodes;

  // initialize cell face counter and cell-node iterator
  unsigned cf_counter = 0;
  auto cn_first = cell_to_node_linkage.begin();

  // convert cell-node linkage to map of cell face to
  for (unsigned cell = 0; cell < num_cells; ++cell) {
    for (unsigned face = 0; face < num_faces_per_cell[cell]; ++face) {

      // convert iterator to node indices to set of node indices
      cface_to_nodes[cf_counter] =
          std::set<unsigned>(cn_first, cn_first + num_nodes_per_face_per_cell[cf_counter]);

      // increment iterator and counter
      cn_first += num_nodes_per_face_per_cell[cf_counter];
      cf_counter++;
    }
  }

  Check(cn_first == cell_to_node_linkage.end());

  // (2) create a map of node-sets to cells

  std::map<std::set<unsigned>, std::vector<unsigned>> nodes_to_cells;

  // reset cf_counter
  cf_counter = 0;

  for (unsigned cell = 0; cell < num_cells; ++cell) {
    for (unsigned face = 0; face < num_faces_per_cell[cell]; ++face) {

      // invert the map
      nodes_to_cells[cface_to_nodes[cf_counter]].push_back(cell);

      // increment counter
      cf_counter++;
    }
  }

  // (3) create maps of nodes to boundary faces (sides) and parallel faces

  std::map<std::set<unsigned>, unsigned> nodes_to_side =
      compute_node_vec_indx_map(side_node_count, side_to_node_linkage);

  std::map<std::set<unsigned>, unsigned> nodes_to_ghost =
      compute_node_vec_indx_map(ghost_cell_type, ghost_cell_to_node_linkage);

  // (4) create cell-to-cell, cell-to-side, cell-to-ghost-cell linkage

  // reset cf_counter and cell-node iterator
  cf_counter = 0;
  cn_first = cell_to_node_linkage.begin();

  for (unsigned cell = 0; cell < num_cells; ++cell) {
    for (unsigned face = 0; face < num_faces_per_cell[cell]; ++face) {

      // initialize this face to not having a condition
      bool has_face_cond = false;

      // get the node set for this cell and face
      const std::set<unsigned> &node_set = cface_to_nodes[cf_counter];

      // get the cells associated with this cell face from the nodes
      const std::vector<unsigned> &cells = nodes_to_cells[node_set];

      Check(cells.size() >= 1);
      Check(cells.size() <= 2);

      // get ordered node vector from cell_to_node_linkage
      const std::vector<unsigned> node_vec(cn_first,
                                           cn_first + num_nodes_per_face_per_cell[cf_counter]);

      // check how many cells are associated with the face
      if (cells.size() == 2) {

        // get neighbor cell index
        const unsigned oth_cell = cell == cells[0] ? cells[1] : cells[0];

        Check(oth_cell != cell);

        // add to cell-cell linkage
        cell_to_cell_linkage[cell].push_back(std::make_pair(oth_cell, node_vec));

        // a neighbor cell was found
        has_face_cond = true;
      }

      // check if a boundary/side exists for this node set
      if (nodes_to_side.find(node_set) != nodes_to_side.end()) {

        // populate cell-boundary face layout
        cell_to_side_linkage[cell].push_back(std::make_pair(nodes_to_side[node_set], node_vec));

        has_face_cond = true;
      }

      // check if a parallel face exists for this node set
      if (nodes_to_ghost.find(node_set) != nodes_to_ghost.end()) {

        // populate cell-parallel face layout
        cell_to_ghost_cell_linkage[cell].push_back(
            std::make_pair(nodes_to_ghost[node_set], node_vec));

        has_face_cond = true;
      }

      // make face a boundary if no face conditions have been found
      if (!has_face_cond) {

        // augment side flags with vacuum b.c.
        side_set_flag.push_back(0);

        // augment side-node count
        m_side_node_count.push_back(num_nodes_per_face_per_cell[cf_counter]);
        Check(m_side_node_count.size() == side_set_flag.size());

        // augment side-node linkage
        m_side_to_node_linkage.insert(m_side_to_node_linkage.begin(), node_vec.begin(),
                                      node_vec.end());

        // augment cell-side linkage
        cell_to_side_linkage[cell].push_back(
            std::make_pair(static_cast<unsigned>(m_side_node_count.size() - 1), node_vec));
      }

      // increment iterator and counter
      cn_first += num_nodes_per_face_per_cell[cf_counter];
      cf_counter++;
    }
  }

  Ensure(cn_first == cell_to_node_linkage.end());
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build a map of node vectors to indices map for boundary layouts.
 *
 * Note: the ordering of the nodes in the mesh constructor must match the node ordering of the
 * corresponding (local) cell face.
 *
 * \param[in] indx_type vector of number of nodes, subscripted by index.
 * \param[in] indx_to_node_linkage serial map of index to node indices.
 * \return a map of node index to vector of indexes adjacent to the node.
 */
std::map<std::set<unsigned>, unsigned>
Draco_Mesh::compute_node_vec_indx_map(const std::vector<unsigned> &indx_type,
                                      const std::vector<unsigned> &indx_to_node_linkage) const {

  // map to return
  std::map<std::set<unsigned>, unsigned> nodes_to_indx_map;

  // generate map
  const size_t num_indxs = indx_type.size();
  auto i2n_first = indx_to_node_linkage.begin();
  for (unsigned indx = 0; indx < num_indxs; ++indx) {

    // extract the node vector
    const std::set<unsigned> node_vec(i2n_first, i2n_first + indx_type[indx]);

    // set the node vector as the key and index as the value
    nodes_to_indx_map[node_vec] = indx;

    // increment iterator
    i2n_first += indx_type[indx];
  }

  Ensure(i2n_first == indx_to_node_linkage.end());

  return nodes_to_indx_map;
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Build a cell-to-cell linkage across corners
 *
 * \param[in] num_faces_per_cell number of faces per cell.
 * \param[in] cell_to_node_linkage serial map of cell to face to node indices.
 * \param[in] num_nodes_per_face_per_cell number of nodes per face per cell.
 * \param[in] ghost_cell_type number of common vertices per ghost cell (sharing a full face).
 * \param[in] ghost_cell_to_node_linkage vertices in common per ghost cell (sharing a full face).
 * \param[in] global_node_number vector indexed by local node with global node index as values.
 */
void Draco_Mesh::compute_node_to_cell_linkage(
    const std::vector<unsigned> &num_faces_per_cell,
    const std::vector<unsigned> &cell_to_node_linkage,
    const std::vector<unsigned> &num_nodes_per_face_per_cell,
    const std::vector<unsigned> &ghost_cell_type,
    const std::vector<unsigned> &ghost_cell_to_node_linkage,
    const std::vector<unsigned> &global_node_number) {

  // (1a) create map of (single) nodes to set of cells

  std::map<unsigned, std::set<unsigned>> node_to_cells;

  // initialize cell face counter and cell-node iterator
  unsigned cf_counter = 0;
  unsigned cn_counter = 0;

  // convert cell-node linkage to map of cell face to
  for (unsigned cell = 0; cell < num_cells; ++cell) {
    for (unsigned face = 0; face < num_faces_per_cell[cell]; ++face) {
      for (unsigned node = 0; node < num_nodes_per_face_per_cell[cf_counter]; ++node) {

        // add cell to the set for the node at this point in the cell-node linkage vector
        node_to_cells[cell_to_node_linkage[cn_counter]].insert({cell});

        // each

        // increment cell node counter
        cn_counter++;
      }

      // increment cell face counter
      cf_counter++;
    }
  }

  Check(cn_counter == safe_convert_from_size_t(cell_to_node_linkage.size()));

  // (1b) convert set of cell per node to vector ...

  for (unsigned node = 0; node < num_nodes; ++node) {

    // get the set of local cells for the node
    const std::set<unsigned> &nc_set = node_to_cells.at(node);

    // convert set to vector
    node_to_cell_linkage[node] = std::vector<unsigned>(nc_set.begin(), nc_set.end());

    // each node must be adjacent to at least one local cell
    Check(node_to_cell_linkage[node].size() >= 1);
  }

  // avoid populating ghost node map if there are no faces that go off rank
  if (ghost_cell_type.size() == 0)
    return;

  //----------------------------------------------------------------------------------------------//
  // When domain-decomposed, the following creates a map of local nodes to all ghost cells.
  // The procedure makes use of existing ghost data across cell faces, as follows:
  //
  // (1) create global node-local cell map with: (local-global node) and (local node-local cell)
  //
  // (2) split map into serialized vectors of local cells and global nodes, over the same index
  //
  // (3) mpi_allgatherv each ranks local cell/global node vectors
  //
  // (4) merged vectors per rank from (3) to per rank vector of map of global node to local cells
  //
  // (5) for a rank (i), compare each other rank (j) ghost global node list
  //
  // (6) if rank (i) has common global ghost nodes with rank (j), add the local cell indices of
  //     rank (j) to each corresponding rank (i) local node index in rank (i)'s dual layout
  //
  //----------------------------------------------------------------------------------------------//

  // create map of global node to vector of adjacent ranks
  std::map<unsigned, std::vector<unsigned>> global_node_to_local_cells;

  // first append this ranks index onto the map
  // initialize cell face counter and cell-node iterator
  unsigned gcn_counter = 0;

  // short-cut to number of ghost faces (i.e. cells across a face on a rank boundary)
  const unsigned num_ghost_cells = safe_convert_from_size_t(ghost_cell_type.size());

  // create the pre-comm map of global ghost nodes to local cells
  for (unsigned ghost = 0; ghost < num_ghost_cells; ++ghost) {
    for (unsigned ghost_node = 0; ghost_node < ghost_cell_type[ghost]; ++ghost_node) {

      // get global index for this node
      const unsigned local_node = ghost_cell_to_node_linkage[gcn_counter];
      const unsigned global_node = global_node_number[local_node];

      // set initial rank and local cell listing at this global node
      global_node_to_local_cells[global_node] = node_to_cell_linkage.at(local_node);

      // increment ghost-cell-node-linkage counter
      gcn_counter++;
    }
  }

  // create a serial array of local cell indices (for cells with at least one off-rank vertex)
  std::vector<unsigned> cells_per_serial;
  for (const auto &global_node_cell_pair : global_node_to_local_cells) {
    cells_per_serial.insert(cells_per_serial.end(), global_node_cell_pair.second.begin(),
                            global_node_cell_pair.second.end());
  }

  // get the serialized vector size
  const size_t num_serial = cells_per_serial.size();

  // initialize a serial array of global node indices
  std::vector<unsigned> global_node_per_serial(num_serial);

  // map global ghost node indices to serial index
  size_t serial_count = 0;
  for (const auto &global_node_cell_pair : global_node_to_local_cells) {

    // get the number of local cells for this node
    const size_t num_node_cells = global_node_cell_pair.second.size();

    // set to the global node at the serial index for this local cell
    for (size_t node_cell = 0; node_cell < num_node_cells; ++node_cell)
      global_node_per_serial[serial_count + node_cell] = global_node_cell_pair.first;

    // increment count over serial index
    serial_count += num_node_cells;
  }

  // check that serial vector has been filled
  Check(serial_count == num_serial);

  // gather the local cell indices per serial index per rank
  std::vector<std::vector<unsigned>> cells_per_serial_per_rank;
  rtt_c4::indeterminate_allgatherv(cells_per_serial, cells_per_serial_per_rank);

  // get the number of ranks in this communicator group
  const unsigned num_ranks = rtt_c4::nodes();

  // sanity-check gatherv communicator group
  Check(cells_per_serial_per_rank.size() == num_ranks);

  // resize gather target for global nodes since the sizes are determined from previous gather
  std::vector<std::vector<unsigned>> global_node_per_serial_per_rank(num_ranks);
  for (unsigned rank = 0; rank < num_ranks; ++rank)
    global_node_per_serial_per_rank[rank].resize(cells_per_serial_per_rank[rank].size());

  // gather the global ghost node indices per serial index per rank
  rtt_c4::determinate_allgatherv(global_node_per_serial, global_node_per_serial_per_rank);

  // merge global_node_per_serial_per_rank and cells_per_serial_per_rank into map per rank
  std::vector<std::map<unsigned, std::vector<unsigned>>> global_node_cell_map_per_rank(num_ranks);
  for (unsigned rank = 0; rank < num_ranks; ++rank) {

    // short-cut to serialized vector size from rank
    const size_t num_serial_on_rank = cells_per_serial_per_rank[rank].size();
    Check(num_serial_on_rank == global_node_per_serial_per_rank[rank].size());

    // generate map for the rank
    for (size_t i = 0; i < num_serial_on_rank; ++i) {

      // short-cut to key (global node index) and value modifier (local_cell on rank)
      const unsigned global_node = global_node_per_serial_per_rank[rank][i];
      const unsigned local_cell = cells_per_serial_per_rank[rank][i];

      // accumulate local cells (for rank) adjacent to this global node
      global_node_cell_map_per_rank[rank][global_node].push_back(local_cell);
    }
  }

  // invert local-to-global node index layout
  std::map<unsigned, unsigned> global_to_local_node;
  for (unsigned node = 0; node < num_nodes; ++node)
    global_to_local_node[global_node_number[node]] = node;

  // get this (my) rank
  const unsigned my_rank = rtt_c4::node();

  // generate dual ghost layout
  for (unsigned rank = 0; rank < num_ranks; ++rank) {

    // exclude this rank
    if (rank == my_rank)
      continue;

    // short-cut to global node vectors to sort and set-intersect
    std::vector<unsigned> &v1 = global_node_per_serial_per_rank[my_rank];
    std::vector<unsigned> &v2 = global_node_per_serial_per_rank[rank];

    // sort the global node vectors
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    // remove duplicates
    auto last_v1 = std::unique(v1.begin(), v1.end());
    v1.erase(last_v1, v1.end());
    auto last_v2 = std::unique(v2.begin(), v2.end());
    v2.erase(last_v2, v2.end());

    // perform set intersection of my rank global nodes with this rank's global nodes
    std::vector<unsigned> rank_my_rank_node_intersect;
    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
                          std::back_inserter(rank_my_rank_node_intersect));

    // map cell and rank from each common global node index to my rank's dual ghost layout
    for (auto gl_node : rank_my_rank_node_intersect) {

      // get the local node index from the global-local-node map
      const unsigned node = global_to_local_node.at(gl_node);

      // append each local-cell-rank pair to dual ghost layout
      for (auto local_cell : global_node_cell_map_per_rank[rank].at(gl_node))
        node_to_ghost_cell_linkage[node].push_back(std::make_pair(local_cell, rank));
    }
  }

  // since this mesh was constructed with ghost data, the resulting map must have non-zero size
  Ensure(node_to_ghost_cell_linkage.size() > 0);
}

} // end namespace rtt_mesh

//------------------------------------------------------------------------------------------------//
// end of mesh/Draco_Mesh.cc
//------------------------------------------------------------------------------------------------//
