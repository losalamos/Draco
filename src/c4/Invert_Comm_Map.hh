//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Invert_Comm_Map.hh
 * \author Mike Buksas
 * \date   Mon Nov 19 10:09:10 2007
 * \brief  
 * \note   Copyright (C) 2007-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_Invert_Comm_Map_hh
#define c4_Invert_Comm_Map_hh

#include "global.hh"
#include <vector>
#include <iterator> // back_inserter

namespace rtt_c4
{

DLL_PUBLIC void master_impl(const std::vector<int> &, std::vector<int>&);
DLL_PUBLIC void slave_impl( const std::vector<int> &, std::vector<int>&);

//---------------------------------------------------------------------------//
/**
 * \brief Invert the contents of a one-to-many mapping between nodes.
 *
 * The argument on each node is a container whose stored values are other
 * nodes. The result is collection of the same type containing the nodes which
 * referred to the local node.
 *
 * E.g. if the argument contains "communicate to" node values. The result
 * contains the "receive from" node values.
 */
template <typename T>
void invert_comm_map(const T& to_values, T& from_values)
{
    const int node = rtt_c4::node();

    // Copy the provided container to a std::vector<int>
    std::vector<int> to_data, from_data;
    to_data.insert(to_data.end(), to_values.begin(), to_values.end());


    if (node == 0)
        master_impl(to_data, from_data);
    else
        slave_impl(to_data, from_data);

    // Append the results to the end of the provided container.
    std::copy(from_data.begin(), from_data.end(), std::back_inserter(from_values));
}

//---------------------------------------------------------------------------//
/**
 * \brief Specialized version of invert_comm_map for std::vector<int> which
 * avoids data copy operations.
 * 
 */
inline void invert_comm_map(const std::vector<int>& to_values,
                            std::vector<int>& from_values)
{
    const int node = rtt_c4::node();

    if (node == 0)
        master_impl(to_values, from_values);
    else
        slave_impl(to_values, from_values);

}

} // end namespace rtt_c4

#endif // c4_Invert_Comm_Map_hh

//---------------------------------------------------------------------------//
// end of c4/Invert_Comm_Map.hh
//---------------------------------------------------------------------------//
