//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   RTT_Format_Reader/RTT_Mesh_Reader.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Mesh_Reader library.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef rtt_RTT_Format_Reader_RTT_Mesh_Reader_hh
#define rtt_RTT_Format_Reader_RTT_Mesh_Reader_hh

#include "RTT_Format_Reader.hh"
#include "meshReaders/Mesh_Reader.hh"

namespace rtt_RTT_Format_Reader {
//============================================================================//
/*!
 * \class RTT_Mesh_Reader
 *
 * \brief An input routine to parse an RTT Format mesh file using the DRACO
 *        meshReaders standard interface.
 *
 *\sa The RTT_Mesh_Reader class is a derived-type of the Mesh_Reader abstract
 *    base class specified in the meshReaders package. Packages using the
 *    RTT_Mesh_Reader should thus include the RTT_Mesh_Reader.hh decorator file
 *    located in the meshReaders directory to resolve the namespace.
 *    RTT_Mesh_Reader contains a RTT_Format_Reader as a private data member, so
 *    none of the RTT_Format_Reader class public accessor functions are
 *    accessible.
 */
//============================================================================//

class DLL_PUBLIC_RTT_Format_Reader RTT_Mesh_Reader
    : public rtt_meshReaders::Mesh_Reader {
  // NESTED CLASSES AND TYPEDEFS
  typedef std::string string;
  typedef std::set<int> set_int;
  typedef std::vector<int> vector_int;
  typedef std::vector<std::vector<int>> vector_vector_int;
  typedef std::vector<std::vector<std::vector<int>>> vector_vector_vector_int;
  typedef std::vector<std::vector<double>> vector_vector_dbl;
  typedef std::set<unsigned> set_uint;
  typedef std::vector<unsigned> vector_uint;
  typedef std::vector<std::vector<unsigned>> vector_vector_uint;
  typedef std::vector<std::vector<std::vector<unsigned>>>
      vector_vector_vector_uint;

  // DATA

private:
  std::shared_ptr<RTT_Format_Reader> rttMesh;
  std::vector<std::shared_ptr<rtt_mesh_element::Element_Definition>>
      element_defs;
  std::vector<rtt_mesh_element::Element_Definition::Element_Type> element_types;
  std::vector<rtt_mesh_element::Element_Definition::Element_Type>
      unique_element_types;

public:
  /*!
   * \brief Constructs an RTT_Mesh_Reader class object.
   * \param RTT_File Mesh file name.
   */
  RTT_Mesh_Reader(const string &RTT_File)
      : rttMesh(new RTT_Format_Reader(RTT_File)), element_defs(),
        element_types(), unique_element_types() {
    transform2CGNS();
  }

  //! Destroys an RTT_Mesh_Reader class object
  ~RTT_Mesh_Reader() {}

  // ACCESSORS

  // Virutal accessor function definitions based on the Mesh_Readers abstract
  // base class.

  /*!
   * \brief Returns the coordinate values for each of the nodes.
   * \return The coordinate values for the nodes.
   */
  virtual vector_vector_dbl get_node_coords() const {
    return rttMesh->get_nodes_coords();
  }

  /*!
   * \brief Returns the problem coordinate units (e.g, cm).
   * \return Coordinate units.
   */
  virtual string get_node_coord_units() const {
    return rttMesh->get_dims_coor_units();
  }

  /*!
   * \brief Returns the topological dimenstion (1, 2 or 3).
   * \return Topological dimension.
   */
  virtual size_t get_dims_ndim() const { return rttMesh->get_dims_ndim(); }

  size_t get_dims_ncells() const { return rttMesh->get_dims_ncells(); }

  size_t get_dims_nsides() const { return rttMesh->get_dims_nsides(); }

  virtual vector_vector_uint get_element_nodes() const;

  /*!
   * \brief Returns the element (i.e., sides and cells) types (e.g., TRI_3 and
   *        TETRA_4).
   * \return Element definitions.
   */
  virtual std::vector<rtt_mesh_element::Element_Definition::Element_Type>
  get_element_types() const {
    return element_types;
  }

  virtual std::vector<std::shared_ptr<rtt_mesh_element::Element_Definition>>
  get_element_defs() const {
    return element_defs;
  }

  /*!
   * \brief Returns the unique element types (e.g., TRI_3 and TETRA_4) that are
   *        defined in the mesh file.
   * \return Element definitions.
   */
  virtual std::vector<rtt_mesh_element::Element_Definition::Element_Type>
  get_unique_element_types() const {
    return unique_element_types;
  }

  virtual std::map<string, set_uint> get_node_sets() const;

  virtual std::map<string, set_uint> get_element_sets() const;
  /*!
   * \brief Returns the mesh file title.
   * \return Mesh file title.
   */
  virtual string get_title() const { return rttMesh->get_header_title(); }

  virtual bool invariant() const;

  // IMPLEMENTATION

private:
  void transform2CGNS(void);
};

} // end namespace rtt_RTT_Format_Reader

#endif // rtt_RTT_Format_Reader_RTT_Mesh_Reader_hh

//----------------------------------------------------------------------------//
// end of RTT_Format_Reader/RTT_Mesh_Reader.hh
//----------------------------------------------------------------------------//
