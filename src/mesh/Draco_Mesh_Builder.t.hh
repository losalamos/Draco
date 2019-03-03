//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Draco_Mesh_Builder.t.hh
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Tuesday, Jul 03, 2018, 11:26 am
 * \brief  Draco_Mesh_Builder class header file.
 * \note   Copyright (C) 2018-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Draco_Mesh.hh"
#include "Draco_Mesh_Builder.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <iostream>

namespace rtt_mesh {

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Draco_Mesh_Builder constructor.
 *
 * \param[in] reader_ shared pointer to a mesh reader object
 */
template <typename FRT>
Draco_Mesh_Builder<FRT>::Draco_Mesh_Builder(std::shared_ptr<FRT> reader_)
    : reader(reader_) {
  Require(reader_ != nullptr);
  // \todo remove constraint of 2 dimensions
  Insist(reader->get_numdim() == 2, "Mesh must be 2D.");
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Build a Draco_Mesh object
 *
 * \param[in] geometry enumeration of mesh geometry
 *
 * \return shared pointer to the Draco_Mesh object
 */
template <typename FRT>
std::shared_ptr<Draco_Mesh>
Draco_Mesh_Builder<FRT>::build_mesh(rtt_mesh_element::Geometry geometry) {

  // \todo: Should geometry be to rtt format parsing?
  // \todo: Generate ghost node and cell data for domain-decomposed meshes.

  Require(geometry != rtt_mesh_element::END_GEOMETRY);
  // \todo: Eventually allow spherical geometry
  Require(geometry != rtt_mesh_element::SPHERICAL);

  // >>> GENERATE MESH CONSTRUCTOR ARGUMENTS

  // get the number of dimensions
  unsigned dimension = reader->get_numdim();
  // \todo: Eventually allow dim = 1, 3
  Check(dimension == 2);

  // get the number of cells
  size_t num_cells = reader->get_numcells();

  // generate the cell type vector
  size_t cn_linkage_size = 0;
  std::vector<unsigned> cell_type(num_cells);
  for (size_t cell = 0; cell < num_cells; ++cell) {

    // for Draco_Mesh, cell_type is number of nodes
    cell_type[cell] = reader->get_celltype(cell);

    // increment size of cell-to-node linkage array
    cn_linkage_size += cell_type[cell];
  }

  // \todo: Can the cell defs past num_cells - 1 be checked as invalid?

  // generate the cell-to-node linkage
  std::vector<unsigned> cell_to_node_linkage;
  cell_to_node_linkage.reserve(cn_linkage_size);
  for (size_t cell = 0; cell < num_cells; ++cell) {

    // insert the vector of node indices
    const std::vector<unsigned> cell_nodes = reader->get_cellnodes(cell);
    cell_to_node_linkage.insert(cell_to_node_linkage.end(), cell_nodes.begin(),
                                cell_nodes.end());
  }

  // get the number of sides
  size_t num_sides = reader->get_numsides();

  // generate the side node count vector
  size_t sn_linkage_size = 0;
  std::vector<unsigned> side_node_count(num_sides);
  std::vector<unsigned> side_set_flag(num_sides);
  for (size_t side = 0; side < num_sides; ++side) {

    // acquire the number of nodes associated with this side def
    Check(reader->get_sidetype(side) < UINT_MAX);
    side_node_count[side] = static_cast<unsigned>(reader->get_sidetype(side));

    // this is not required in rtt meshes, but is so in Draco_Mesh
    Check(dimension == 2 ? side_node_count[side] == 2 : true);

    // get the 1st side flag associated with this side
    // \todo: What happens when side has no flags?
    side_set_flag[side] = reader->get_sideflag(side);

    // increment size of cell-to-node linkage array
    sn_linkage_size += side_node_count[side];
  }

  // \todo: Can the side defs past num_sides - 1 be checked as invalid?

  // generate the side-to-node linkage
  std::vector<unsigned> side_to_node_linkage;
  side_to_node_linkage.reserve(sn_linkage_size);
  for (size_t side = 0; side < num_sides; ++side) {

    // insert the vector of node indices
    const std::vector<unsigned> side_nodes = reader->get_sidenodes(side);
    side_to_node_linkage.insert(side_to_node_linkage.end(), side_nodes.begin(),
                                side_nodes.end());
  }

  // get the number of nodes
  size_t num_nodes = reader->get_numnodes();

  Check(num_nodes >= num_cells);
  Check(num_nodes <= cn_linkage_size);

  // \todo: add global node numbers to rtt mesh reader?
  // \todo: or, remove global_node_number from mesh constructor?
  // assume domain is not decomposed, for now

  // generate the global node number serialized vector of coordinates
  std::vector<unsigned> global_node_number(num_nodes);
  std::vector<double> coordinates(dimension * num_nodes);
  for (size_t node = 0; node < num_nodes; ++node) {

    // set the "global" node indices
    Check(node < UINT_MAX);
    global_node_number[node] = static_cast<unsigned>(node);

    // get coordinates for this node
    const std::vector<double> node_coord = reader->get_nodecoord(node);

    // populate coordinate vector
    for (unsigned d = 0; d < dimension; ++d)
      coordinates[dimension * node + d] = node_coord[d];
  }

  Remember(auto cn_minmax = std::minmax_element(cell_to_node_linkage.begin(),
                                                cell_to_node_linkage.end()));
  Remember(auto sn_minmax = std::minmax_element(side_to_node_linkage.begin(),
                                                side_to_node_linkage.end()));
  Ensure(*cn_minmax.first >= 0);
  Ensure(*cn_minmax.second < num_nodes);
  Ensure(side_to_node_linkage.size() > 0 ? *sn_minmax.first >= 0 : true);
  Ensure(side_to_node_linkage.size() > 0 ? *sn_minmax.second < num_nodes
                                         : true);

  // >>> CONSTRUCT THE MESH

  std::shared_ptr<Draco_Mesh> mesh(new Draco_Mesh(
      dimension, geometry, cell_type, cell_to_node_linkage, side_set_flag,
      side_node_count, side_to_node_linkage, coordinates, global_node_number));

  return mesh;
}

} // end namespace rtt_mesh

//---------------------------------------------------------------------------//
// end of mesh/Draco_Mesh_Builder.t.hh
//---------------------------------------------------------------------------//
