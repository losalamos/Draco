//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/TET_test_1.hh
 * \author H. Grady Hughes
 * \date   Tue Feb  8 18:13:29 MST 2000
 * \brief  Two-tet pyramid model class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_test_TET_test_1_hh__
#define __mc_test_TET_test_1_hh__

#include <string>
#include <vector>
#include <map>
#include <set>
#include "meshReaders/Element_Definition.hh"

namespace rtt_mc_test
{

using rtt_meshReaders::Element_Definition;

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

bool fail(int line, char *file)
{
    std::cout << "Test: failed on line " << line << " in " << file
          << std::endl;
    return false;
}

//___________________________________________________________________________//
/*!
 * \class TET_test_1
 * \brief Two-tet pyramid model for testing TET_Builder.
 *
 * TET_test_1 contains no private data.  The specific test case is hard-wired,
 * and delivered to TET_Builder by the public interface functions.
 */
// revision history:
// -----------------
//  0)   Original: Committed 2000-02-09
//  1) 2000-04-10: Rewritten to be consistent with new meshReader classes.
//  2) 2000-04-26: Now uses get_element_types() from meshReader class.
//  3) 2000-05-03: TET_Builder, TET_Mesh, and their test files now use the
//                 get_node_coord_units(), get_node_sets(), get_element_sets(),
//                 and get_title() services of the Mesh_Reader base class.
//                 At the top level (TET_Mesh), the get_element_sets() services
//                 will later be replaced by side- and cell-specific data
//                 structures.
//  4) 2000-06-08: Information from the interface service get_element_sets()
//                 is now converted to two separate maps, side_sets and
//                 cell_sets, and used to initialize data members of the
//                 TET_Mesh class.  The TET_Mesh class no longer has knowledge
//                 of element_sets.  New diagnostic functions print_node_sets,
//                 print_side_sets, and print_cell_sets are added to TET_Mesh.
//
//___________________________________________________________________________//

class TET_test_1
{

 public:

    //! Typedef for a standard set of integers.
    typedef std::set<int> SetInt;

    //! Typedef for a map linking strings to sets of integers.
    typedef std::map< std::string, SetInt > MAP_String_SetInt;

    // Default constructor.  No need to initialize anything.
    TET_test_1() { }

    // TET_Builder no longer requires this function.
    // std::string get_coord_system();

    std::string get_title()
    {
        return std::string("Test 2-tet pyramid in RTT-format mesh.");
    }

    std::vector< std::vector<double> > get_node_coords()
    {
        std::vector< std::vector<double> > node_coor;
        node_coor.resize(5);

        node_coor[0].push_back(0.0);
        node_coor[0].push_back(0.0);
        node_coor[0].push_back(0.0);

        node_coor[1].push_back(1.0);
        node_coor[1].push_back(0.0);
        node_coor[1].push_back(0.0);

        node_coor[2].push_back(0.0);
        node_coor[2].push_back(1.0);
        node_coor[2].push_back(0.0);

        node_coor[3].push_back(1.0);
        node_coor[3].push_back(1.0);
        node_coor[3].push_back(0.0);

        node_coor[4].push_back(0.5);
        node_coor[4].push_back(0.5);
        node_coor[4].push_back(1.0);


        return node_coor;
    }

    // TET_Builder no longer requires this function.
    // std::vector<int> get_parent();

    std::vector< std::vector<int> > get_element_nodes()
    {
        std::vector< std::vector<int> > elements_nodes;
        elements_nodes.resize(8);

        elements_nodes[0].push_back(0);
        elements_nodes[0].push_back(4);
        elements_nodes[0].push_back(2);

        elements_nodes[1].push_back(0);
        elements_nodes[1].push_back(1);
        elements_nodes[1].push_back(4);

        elements_nodes[2].push_back(0);
        elements_nodes[2].push_back(2);
        elements_nodes[2].push_back(1);

        elements_nodes[3].push_back(1);
        elements_nodes[3].push_back(2);
        elements_nodes[3].push_back(3);

        elements_nodes[4].push_back(1);
        elements_nodes[4].push_back(3);
        elements_nodes[4].push_back(4);

        elements_nodes[5].push_back(2);
        elements_nodes[5].push_back(4);
        elements_nodes[5].push_back(3);

        elements_nodes[6].push_back(0);
        elements_nodes[6].push_back(1);
        elements_nodes[6].push_back(2);
        elements_nodes[6].push_back(4);

        elements_nodes[7].push_back(2);
        elements_nodes[7].push_back(1);
        elements_nodes[7].push_back(4);
        elements_nodes[7].push_back(3);

        return elements_nodes;
    }

    std::vector<Element_Definition::Element_Type> get_element_types()
    {
        std::vector<Element_Definition::Element_Type> element_types;

        element_types.push_back(Element_Definition::TRI_3);
        element_types.push_back(Element_Definition::TRI_3);
        element_types.push_back(Element_Definition::TRI_3);
        element_types.push_back(Element_Definition::TRI_3);
        element_types.push_back(Element_Definition::TRI_3);
        element_types.push_back(Element_Definition::TRI_3);
        element_types.push_back(Element_Definition::TETRA_4);
        element_types.push_back(Element_Definition::TETRA_4);

        return element_types;
    }

    std::string get_node_coord_units()
    {
        return std::string("cm");
    }

    MAP_String_SetInt get_node_sets()
    {
        MAP_String_SetInt node_sets;

        SetInt s1;
        s1.insert(1);
        s1.insert(2);
        node_sets[std::string("node_type/dudded")] = s1;

        SetInt s2;
        s2.insert(0);
        node_sets[std::string("node_type/interior")] = s2;

        SetInt s3;
        s3.insert(3);
        s3.insert(4);
        node_sets[std::string("node_type/parent")] = s3;

        SetInt s4;
        s4.insert(0);
        s4.insert(1);
        node_sets[std::string("radiation_boundary/reflective")] = s4;

        SetInt s5;
        s5.insert(2);
        s5.insert(3);
        s5.insert(4);
        node_sets[std::string("radiation_boundary/vacuum")] = s5;

        SetInt s6;
        s6.insert(0);
        s6.insert(1);
        s6.insert(2);
        node_sets[std::string("radiation_source/no_source")] = s6;

        SetInt s7;
        s7.insert(3);
        s7.insert(4);
        node_sets[std::string("radiation_source/rad_source")] = s7;

        return node_sets;
    }

    MAP_String_SetInt get_element_sets()
    {
        MAP_String_SetInt element_sets;

        SetInt s8;
        s8.insert(0);
        s8.insert(1);
        s8.insert(2);
        s8.insert(3);
        s8.insert(4);
        s8.insert(5);
        element_sets[std::string("boundary_conditions/reflective")] = s8;

        SetInt s9;
        element_sets[std::string("boundary_conditions/vacuum")] = s9;

        SetInt s10;
        element_sets[std::string("ion_source/src_name1")] = s10;

        SetInt s11;
        s11.insert(6);
        s11.insert(7);
        element_sets[std::string("ion_source/src_name2")] = s11;

        SetInt s12;
        s12.insert(6);
        s12.insert(7);
        element_sets[std::string("material_region/control_rod")] = s12;

        SetInt s13;
        element_sets[std::string("material_region/shield")] = s13;

        return element_sets;
    }

    // TET_Builder no longer requires this function.
    // bool get_submesh();

};  // end class TET_test_1

}   // end namespace rtt_mc_test

#endif  // __mc_test_TET_test_1_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/TET_test_1.hh
//---------------------------------------------------------------------------//
