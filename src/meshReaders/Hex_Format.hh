//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/Hex_Format.hh
 * \author John McGhee
 * \date   Tue Mar  7 08:38:04 2000
 * \brief  Header file for CIC-19 Hex format mesh reader.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_Hex_Format_hh__
#define __meshReaders_Hex_Format_hh__

#include <vector>
#include <string>
#include "Mesh_Reader.hh"

namespace rtt_meshReaders
{
 
//===========================================================================//
/*!
 * \class Hex_Format
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Hex_Format : public rtt_meshReaders::Mesh_Reader
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    std::string meshfile_name;
    static std::string keyword()
    {
	return "cic19_hex_mesh";
    }
    std::string version;
    int npoints;
    int ncells;
    int nvrtx;
    int nvrpf;
    int ndim;
    int nvb_faces;
    int nrb_faces;
    int nmat;
    std::vector<std::vector<double> > point_coords;
    std::vector<std::vector<int> > ipar;
    std::vector<int> imat_index;
    std::vector<int> irgn_vb_index;
    std::vector<std::vector<int> > ipar_vb;
    std::vector<std::vector<int> > ipar_rb;

    std::map<std::string, std::set<int> > node_sets;
    
  public:

    // CREATORS
    
    Hex_Format(std::string filename);

    // Defaulted Hex_Format(const Hex_Format &rhs);
    // Defaulted ~Hex_Format();

    // MANIPULATORS
    
    // Defaulted Hex_Format& operator=(const Hex_Format &rhs);

    // ACCESSORS

    std::vector<std::vector<double> > get_node_coords() const
    {
	return point_coords;
    }

    std::string get_node_coord_units() const
    {
	return "unknown";
    }

     std::map<std::string, std::set<int> > get_node_sets() const
    {
	return node_sets;
    }

     std::string get_title() const
    {
	return "Untitled -- CIC-19 Hex Mesh";
    }

    std::vector<std::vector<int> > get_element_nodes() const;
    bool invariant() const;
    std::map<std::string, std::set<int> > get_element_sets() const;
    std::vector<Element_Definition::Element_Type> get_element_types() const;
  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_meshReaders

#endif                          // __meshReaders_Hex_Format_hh__

//---------------------------------------------------------------------------//
//                              end of meshReaders/Hex_Format.hh
//---------------------------------------------------------------------------//
