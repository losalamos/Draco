//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestHexMeshReader.hh
 * \author John McGhee
 * \date   Thu Mar  9 08:54:59 2000
 * \brief  Header file for the Hex_Mesh_Reader class unit test.
 * \note   Copyright (C) 2002-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_test_TestHexMeshReader_hh__
#define __meshReaders_test_TestHexMeshReader_hh__

#include <string>
#include <map>
#include <set>

namespace rtt_meshReaders
{
class Hex_Mesh_Reader;
}

namespace rtt_meshReaders_test
{

bool check_mesh(const rtt_meshReaders::Hex_Mesh_Reader &mesh, 
                const std::string &testid);
bool check_nodes(const rtt_meshReaders::Hex_Mesh_Reader &mesh, 
                 const std::string &testid);
bool check_node_units(const rtt_meshReaders::Hex_Mesh_Reader &mesh); 
bool check_node_sets(const rtt_meshReaders::Hex_Mesh_Reader &mesh, 
                     const std::string &testid); 
bool check_title(const rtt_meshReaders::Hex_Mesh_Reader &mesh); 
bool check_element_nodes(const rtt_meshReaders::Hex_Mesh_Reader &mesh, 
                         const std::string &testid);
bool check_invariant(const rtt_meshReaders::Hex_Mesh_Reader &mesh);
bool check_element_sets(const rtt_meshReaders::Hex_Mesh_Reader &mesh, 
                        const std::string &testid);
bool check_element_types(const rtt_meshReaders::Hex_Mesh_Reader &mesh, 
                         const std::string &testid);
bool check_unique_element_types(const rtt_meshReaders::Hex_Mesh_Reader 
                                &mesh, const std::string &testid);
bool compare_double(const double &lhs, const double &rhs);
bool check_map(const std::map<std::string, std::set<int> >
               &elmsets, const std::string &name, const int &begin, 
               const int &end);
bool check_get_dims_ndim(const rtt_meshReaders::Hex_Mesh_Reader &mesh, 
                         const std::string &testid);


} // end namespace rtt_meshReaders_test

#endif // __meshReaders_test_TestHexMeshReader_hh__

//---------------------------------------------------------------------------//
// end of meshReaders/test/TestHexMeshReader.hh
//---------------------------------------------------------------------------//
