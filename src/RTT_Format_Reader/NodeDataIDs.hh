//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   RTT_Format_Reader/NodeDataIDs.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/NodeDataIDs class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_NodeDataIDs_hh__
#define __RTT_Format_Reader_NodeDataIDs_hh__

#include "Dims.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the 
 *        node data ids block of the mesh file.
 */
class NodeDataIDs {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;
  typedef std::vector<string> vector_str;

  const Dims &dims;
  vector_str names;
  vector_str units;

public:
  NodeDataIDs(const Dims &dims_)
      : dims(dims_), names(dims.get_nnode_data()),
        units(dims.get_nnode_data()) {}
  ~NodeDataIDs() {}

  void readDataIDs(ifstream &meshfile);

private:
  void readKeyword(ifstream &meshfile);
  void readData(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
   * \brief Returns the specified node_data_id name.
   * \param id_numb node_data_id index number.
   * \return The node_data_id name.
   */
  string get_data_id_name(size_t id_numb) const {
    Insist(id_numb <= dims.get_nnode_data() - 1,
           "Invalid node data id number!");
    return names[id_numb];
  }

  /*!
   * \brief Returns the units associated with the specified node_data_id.
   * \param id_numb node_data_id index number.
   * \return The node_data_id units.
   */
  string get_data_id_units(size_t id_numb) const {
    Insist(id_numb <= dims.get_nnode_data() - 1,
           "Invalid node data id number!");
    return units[id_numb];
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // __RTT_Format_Reader_NodeDataIDs_hh__

//---------------------------------------------------------------------------//
// end of RTT_Format_Reader/NodeDataIDs.hh
//---------------------------------------------------------------------------//
