//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/Mesh_Reader.hh
 * \author John McGhee
 * \date   Fri Feb 25 08:14:54 2000
 * \brief  Header file for the RTT Mesh_Reader base class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_Mesh_Reader_hh__
#define __meshReaders_Mesh_Reader_hh__

#include <vector>
#include <set>
#include <string>
#include <map>
#include "Element_Definition.hh"

namespace rtt_meshReaders
{
 
//===========================================================================//
/*!
 * \class Mesh_Reader
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Mesh_Reader 
{

    // NESTED CLASSES AND TYPEDEFS


    // DATA
    
  public:

    // CREATORS
    
    //Defaulted: Mesh_Reader();
    //Defaulted: Mesh_Reader(const Mesh_Reader &rhs);

    virtual ~Mesh_Reader() 
    {
	//Empty
    }

    // MANIPULATORS
    
    //Defaulted: Mesh_Reader& operator=(const Mesh_Reader &rhs);

    // ACCESSORS

    virtual std::vector<std::vector<double> > get_node_coords() const = 0;
    virtual std::string get_node_coord_units() const = 0;
    virtual std::vector<std::vector<int> > get_element_nodes() const = 0;
    virtual std::vector<Element_Definition::Element_Type> get_element_types() 
        const = 0;
    virtual std::map<std::string, std::set<int> > get_node_sets() const = 0;
    virtual std::map<std::string, std::set<int> > get_element_sets() const = 0;
    virtual std::string get_title() const = 0;
    virtual bool invariant() const = 0;

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_meshReaders

#endif                          // __meshReaders_Mesh_Reader_hh__

//---------------------------------------------------------------------------//
//                              end of meshReaders/Mesh_Reader.hh
//---------------------------------------------------------------------------//
