//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Invert_Comm_Map.cc
 * \author Mike Buksas
 * \date   Mon Nov 19 10:09:11 2007
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "Invert_Comm_Map.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <map>

namespace rtt_c4
{

const int SIZE_CHANNEL = 325;
const int MAP_CHANNEL  = 326;

const int HOST = 0;

typedef std::map<int, std::vector<int> > map_type;

void master_impl(const std::vector<int> &in, std::vector<int> &out)
{

    const int nodes = rtt_c4::nodes();

    // Create a complete map and copy the local map into it.
    map_type in_maps;
    in_maps[HOST] = in;

    // Receive the map data from all other nodes.
    for (int node=1; node<nodes; ++node)
    {
        int size;
        receive(&size, 1, node, SIZE_CHANNEL);
        Check (size >= 0);

        if (size > 0)
        {
            in_maps[node].resize(size);
            receive(&(in_maps[node])[0], size, node, MAP_CHANNEL);
        }

    }


    // Build an out-going map from the contents of the in-going map.
    map_type out_maps;

    for (map_type::const_iterator from_node = in_maps.begin();
         from_node != in_maps.end();
         ++from_node)
    {
        // Grab the from-node value and a reference to the list of to-nodes.
        const int from_node_val          = from_node->first;
        const std::vector<int> &to_nodes = from_node->second;
        
        // For each to-node, add the from-node to it's data.
        for (std::vector<int>::const_iterator to_node = to_nodes.begin();
             to_node != to_nodes.end();
             ++to_node)
        {

            Check (*to_node >= 0);
            Check (*to_node < nodes);

            out_maps[*to_node].push_back(from_node_val);
        }
    }



    // Communicate the results for remote nodes
    for (map_type::const_iterator to_node = out_maps.begin();
         to_node != out_maps.end();
         ++to_node)
    {
        // Grab the to-node value and a reference to the values.
        const int to_node_val              = to_node->first;
        const std::vector<int> &from_nodes = to_node->second;

        if (to_node_val == HOST) continue; // Don't communicate to host.
        
        // Send the size, and the data if size > 0
        const int size = from_nodes.size();
        send(&size, 1, to_node_val, SIZE_CHANNEL);
        if (size > 0)
        {
            send(&from_nodes[0], size, to_node_val, MAP_CHANNEL);
        }

    }

    // Copy the results for the host node.
    std::copy(out_maps[HOST].begin(), out_maps[HOST].end(), back_inserter(out));

}


void slave_impl(const std::vector<int> &in, std::vector<int> &out)
{

    // Send size of in map and, if > 0, contents
    const int size = in.size();
    send(&size, 1, HOST, SIZE_CHANNEL);

    if (size > 0)
    {
        send(&in[0], size, HOST, MAP_CHANNEL);
    }

    
    // Receive results from the host.
    int recv_size;
    receive(&recv_size, 1, HOST, SIZE_CHANNEL);

    if (recv_size > 0)
    {
        out.resize(recv_size);
        receive(&out[0], recv_size, HOST, MAP_CHANNEL);

    }
    
}



} // end namespace rtt_c4

//---------------------------------------------------------------------------//
//                 end of Invert_Comm_Map.cc
//---------------------------------------------------------------------------//
