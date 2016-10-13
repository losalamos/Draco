//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Invert_Comm_Map.cc
 * \author Mike Buksas
 * \date   Mon Nov 19 10:09:11 2007
 * \brief  
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Invert_Comm_Map.hh"
#include "ds++/Assert.hh"
#include <map>

namespace rtt_c4 {

const int SIZE_CHANNEL = 325;
const int MAP_CHANNEL = 326;
const int HOST = 0;

typedef std::map<int, std::vector<int>> map_type;

void icm_master_impl(std::vector<int> const &in, std::vector<int> &out) {
  int const nodes = rtt_c4::nodes();

  // Create a complete map and copy the local map into it.
  map_type in_maps;
  in_maps[HOST] = in;

  // Receive the map data from all other nodes.
  for (int node = 1; node < nodes; ++node) {
    int size = -1;
    receive(&size, 1, node, SIZE_CHANNEL);
    Check(size >= 0);

    if (size > 0) {
      in_maps[node].resize(size);
      receive(&(in_maps[node])[0], size, node, MAP_CHANNEL);
    }
  }

  // Build an out-going map from the contents of the in-going map.
  map_type out_maps;

  for (map_type::const_iterator from_node = in_maps.begin();
       from_node != in_maps.end(); ++from_node) {
    // Grab the from-node value and a reference to the list of to-nodes.
    const int from_node_val = from_node->first;
    const std::vector<int> &to_nodes = from_node->second;

    // For each to-node, add the from-node to it's data.
    for (std::vector<int>::const_iterator to_node = to_nodes.begin();
         to_node != to_nodes.end(); ++to_node) {

      Check(*to_node >= 0);
      Check(*to_node < nodes);

      out_maps[*to_node].push_back(from_node_val);
    }
  }

  // Communicate the results for all remote nodes
  for (int node = 1; node < nodes; ++node) {
    if (out_maps.count(node) && out_maps[node].size() > 0) {
      // Send size and data
      const int size = out_maps[node].size();
      send(&size, 1, node, SIZE_CHANNEL);
      send(&(out_maps[node])[0], size, node, MAP_CHANNEL);
    } else {
      // Send zero for the size and no data.
      const int size = 0;
      send(&size, 1, node, SIZE_CHANNEL);
    }
  }

  // Copy the results for the host node.
  std::copy(out_maps[HOST].begin(), out_maps[HOST].end(), back_inserter(out));
  return;
}

//----------------------------------------------------------------------------//
void icm_slave_impl(std::vector<int> const &in, std::vector<int> &out) {
  // Send size of in map and, if > 0, contents
  int const size = in.size();
  send(&size, 1, HOST, SIZE_CHANNEL);

  if (size > 0) {
    send(&in[0], size, HOST, MAP_CHANNEL);
  }

  // Receive results from the host.
  int recv_size = -1;
  receive(&recv_size, 1, HOST, SIZE_CHANNEL);

  if (recv_size > 0) {
    out.resize(recv_size);
    receive(&out[0], recv_size, HOST, MAP_CHANNEL);
  }
  return;
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// end of Invert_Comm_Map.cc
//---------------------------------------------------------------------------//
