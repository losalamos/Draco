//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/CellDefs.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for CellDefs library.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_CellDefs_hh__
#define __RTT_Format_Reader_CellDefs_hh__

#include "Dims.hh"
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace rtt_RTT_Format_Reader {
// forward declaration.
class CellDefs;

//---------------------------------------------------------------------------//
/*!
 * \brief Controls parsing, storing, and accessing the data contained in the
 *        cell definition block of the mesh file.
 */
//---------------------------------------------------------------------------//
class CellDef {
  /* TYPEDEFS */

  typedef std::ifstream ifstream;
  typedef std::string string;
  typedef std::set<int> set_int;
  typedef std::vector<int> vector_int;
  typedef std::vector<std::vector<int>> vector_vector_int;

  /* DATA */

  const CellDefs &cellDefs;
  const string name;
  size_t nnodes;
  size_t nsides;
  vector_int side_types;
  std::vector<std::vector<size_t>> sides;
  // Add the capability to maintain the sense of the outward normals.
  std::vector<std::vector<size_t>> ordered_sides;
  // Mapping between the old and new cell definition nodes.
  vector_int node_map;

public:
  CellDef(const CellDefs &cellDefs_, const string &name_)
      : cellDefs(cellDefs_), name(name_), nnodes(0), nsides(0),
        side_types(std::vector<int>()),
        sides(std::vector<std::vector<size_t>>()), ordered_sides(0),
        node_map() { /* empty */
  }

  ~CellDef(void) { /*empty*/
  }

  void readDef(ifstream &meshfile);
  void
  redefineCellDef(vector_int const &new_side_types_,
                  std::vector<std::vector<size_t>> const &new_ordered_sides,
                  int const ndim);

public:
  /*!
 * \brief  Returns the cell definition name.
 * \return The cell definition name.
 */
  string get_name() const { return name; }
  /*!
 * \brief  Returns the number of nodes associated with the cell definition.
 * \return The number of nodes comprising the cell definition.
 */
  size_t get_nnodes() const { return nnodes; }
  /*!
 * \brief  Returns the number of sides associated with the cell definition.
 * \return The number of sides comprising the cell definition.
 */
  size_t get_nsides() const { return nsides; }
  /*!
 * \brief   Returns the side type number associated with the cell definition
 *          specified side index.
 * \param s Side index number.
 * \return  The side type number.
 */
  int get_side_types(size_t s) const { return side_types[s]; }

  vector_int get_all_side_types() const { return side_types; }
  std::vector<std::vector<size_t>> get_all_sides() const { return sides; }
  std::vector<std::vector<size_t>> get_all_ordered_sides() const {
    return ordered_sides;
  }

  /*!
     * \brief Returns the side definition of the specified side index of this
     *        cell definition with the returned cell-node indexes in sorted
     *        order.
     * \param s Side index number.
     * \return The side definition (i.e., the cell-node indexes that comprise
     *        the side).
     */
  std::vector<size_t> const &get_side(size_t s) const { return sides[s]; }
  /*!
     * \brief Returns the side definition of the specified side index of this
     * cell
     *        definition with the returned cell-node indexes ordered to preserve
     *        the right hand rule for the outward-directed normal.
     * \param s Side index number.
     * \return The side definition (i.e., the cell-node indexes that comprise
     * the
     *         side).
     */
  std::vector<size_t> const &get_ordered_side(size_t s) const {
    return ordered_sides[s];
  }

  //---------------------------------------------------------------------------//

  /*!
     * \brief Returns the new nodes map when cell redefinition has been
     *        performed.
     * \return New nodes map.
     */
  const vector_int &get_node_map() const { return node_map; }
};

//---------------------------------------------------------------------------//
/*!
 * \brief Controls parsing, storing, and accessing the data contained in the
 *        cell definition block of the mesh file.
 */
class CellDefs {
  // typedefs
  typedef std::ifstream ifstream;
  typedef std::string string;
  typedef std::set<int> set_int;
  typedef std::vector<int> vector_int;
  typedef std::vector<std::vector<int>> vector_vector_int;
  typedef std::vector<std::vector<std::vector<int>>> vector_vector_vector_int;

  const Dims &dims;
  std::vector<std::shared_ptr<CellDef>> defs;
  bool redefined;

public:
  CellDefs(const Dims &dims_)
      : dims(dims_), defs(dims.get_ncell_defs()), redefined(false) { /*empty*/
  }

  ~CellDefs() { /*empty*/
  }

  void readCellDefs(ifstream &meshfile);
  void redefineCellDefs(
      vector_vector_int const &cell_side_types,
      std::vector<std::vector<std::vector<size_t>>> const &cell_ordered_sides);

private:
  void readKeyword(ifstream &meshfile);
  void readDefs(ifstream &meshfile);
  void readEndKeyword(ifstream &meshfile);

public:
  /*!
 * \brief Returns the name of the specified cell definition.
 * \param i Cell definition index number.
 * \return The cell definition name.
 */
  string get_name(unsigned i) const {
    Require(i < defs.size());
    return defs[i]->get_name();
  }
  /*!
 * \brief Returns the specified cell definition.
 * \param i Cell definition index number.
 * \return The cell definition.
 */
  const CellDef &get_cell_def(size_t i) const { return *(defs[i]); }
  std::shared_ptr<CellDef> get_def(size_t i) const { return defs[i]; }
  /*!
 * \brief Returns the number of nodes associated with the specified cell
 *        definition.
 * \param i Cell definition index number.
 * \return The number of nodes comprising the cell definition.
 */
  size_t get_nnodes(size_t i) const { return defs[i]->get_nnodes(); }
  /*!
 * \brief Returns the number of sides associated with the specified cell
 *        definition.
 * \param i Cell definition index number.
 * \return The number of sides comprising the cell definition.
 */
  size_t get_nsides(size_t i) const { return defs[i]->get_nsides(); }
  /*!
 * \brief Returns the side type number associated with the specified side
 *        index and cell definition.
 * \param i Cell definition index number.
 * \param s Side index number.
 * \return The side type number.
 */
  int get_side_types(size_t i, size_t s) const {
    return defs[i]->get_side_types(s);
  }
  /*!
     * \brief Returns the side definition associated with the specified cell
     *        definition and side index with the returned cell-node indexes in
     *        sorted order.
     * \param i Cell definition index number.
     * \param s Side index number.
     * \return The side definition (i.e., the cell-node indexes that comprise
     *        the side).
     */
  std::vector<size_t> const &get_side(int i, int s) const {
    return defs[i]->get_side(s);
  }
  /*!
     * \brief Returns the side definition associated with the specified cell
     *        definition and side index with the returned cell-node indexes
     *        ordered to preserve the right hand rule for the outward-directed
     *        normal.
     * \param i Cell definition index number.
     * \param s Side index number.
     * \return The side definition (i.e., the cell-node indexes that comprise
     *        the side).
     */
  std::vector<size_t> const &get_ordered_side(int i, int s) const {
    return defs[i]->get_ordered_side(s);
  }
  /*!
 * \brief Returns the number of unique cell type definitions.
 * \return The number of cell definitions.
 */
  //    int get_ncell_defs() const { return dims.get_ncell_defs(); }
  /*!
     * \brief Returns the status of the flag indicating that the cell
     * definitions
 *        have been redefined.
 * \return The status of the redefined flag.
 */
  bool get_redefined() const { return redefined; }
  /*!
 * \brief Returns the new node map for the specified cell definition when
 *        redefinition has been performed.
 * \param cell_def Cell definition index.
 * \return New cell definition node map.
 */
  const vector_int &get_node_map(int cell_def) const {
    return defs[cell_def]->get_node_map();
  }
};

} // end namespace rtt_RTT_Format_Reader

#endif // __RTT_Format_Reader_CellDefs_hh__

//---------------------------------------------------------------------------//
// end of RTT_Format_Reader/CellDefs.hh
//---------------------------------------------------------------------------//
