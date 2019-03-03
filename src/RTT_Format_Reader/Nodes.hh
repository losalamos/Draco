//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   RTT_Format_Reader/Nodes.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Nodes class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_Nodes_hh__
#define __RTT_Format_Reader_Nodes_hh__

#include "NodeFlags.hh"

namespace rtt_RTT_Format_Reader {
/*!
 * \brief Controls parsing, storing, and accessing the data specific to the 
 *        nodes block of the mesh file.
 */
class Nodes {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;
  typedef std::vector<int> vector_int;
  typedef std::vector<std::vector<int>> vector_vector_int;
  typedef std::vector<double> vector_dbl;
  typedef std::vector<std::vector<double>> vector_vector_dbl;

  const NodeFlags &nodeFlags;
  const Dims &dims;
  vector_vector_dbl coords;
  vector_int parents;
  vector_vector_int flags;

public:
  Nodes(const NodeFlags &nodeFlags_, const Dims &dims_)
      : nodeFlags(nodeFlags_), dims(dims_),
        coords(dims.get_nnodes(), vector_dbl(dims.get_ndim())),
        parents(dims.get_nnodes()),
        flags(dims.get_nnodes(), vector_int(dims.get_nnode_flag_types())) {}
  ~Nodes() {}

  void readNodes(ifstream &meshfile);
  static int readNextInt(ifstream &meshfile);

private:
  void readKeyword(ifstream &meshfile);
  void readData(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
   * \brief Returns the coordinate values for each of the nodes.
   * \return The coordinate values for the nodes.
   */
  vector_vector_dbl get_coords() const { return coords; }

  /*!
   * \brief Returns all of the coordinate values for the specified node.
   * \param node_numb Node number.
   * \return The node coordinate values.
   */
  vector_dbl get_coords(size_t node_numb) const { return coords[node_numb]; }

  /*!
   * \brief Returns the coordinate value for the specified node and direction 
   *        (i.e., x, y, and z).
   * \param node_numb Node number.
   * \param coord_index Coordinate index number (x = 0, y = 1, z = 2).
   * \return The node coordinate value.
   */
  double get_coords(size_t node_numb, size_t coord_index) const {
    return coords[node_numb][coord_index];
  }

  /*!
   * \brief Returns the node parent for the specified node.
   * \param node_numb Node number.
   * \return The node parent.
   */
  int get_parents(size_t node_numb) const { return parents[node_numb]; }

  /*!
   * \brief Returns the node flag for the specified node and flag index.
   * \param node_numb Node number.
   * \param flag_numb Node flag index.
   * \return The node flag.
   */
  int get_flags(size_t node_numb, size_t flag_numb) const {
    return flags[node_numb][flag_numb];
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // __RTT_Format_Reader_Nodes_hh__

//---------------------------------------------------------------------------//
// end of RTT_Format_Reader/Nodes.hh
//---------------------------------------------------------------------------//
