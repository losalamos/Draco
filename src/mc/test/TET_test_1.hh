//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/TET_test_1.hh
 * \author H. Grady Hughes
 * \date   Tue Feb  8 18:13:29 MST 2000
 * \brief  Two-tet pyramid model class header file.
 */
//---------------------------------------------------------------------------//

#ifndef __mc_test_TET_test_1_hh__
#define __mc_test_TET_test_1_hh__

#include <string>
#include <vector>
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
//  0) original : Committed 2000-02-09
//  1) 2000-04-10 : Rewritten to be consistent with new meshReader classes.
//___________________________________________________________________________//

class TET_test_1
{

 public:

    // Default constructor.  No need to initialize anything.
    TET_test_1() { }

    // TET_Builder no longer requires this function.
    // std::string get_coord_system();

    std::vector< std::vector<double> > get_node_coords()
    {
        std::vector< std::vector<double> > nodes_coor;
        nodes_coor.resize(5);

        nodes_coor[0].push_back(0.0);
        nodes_coor[0].push_back(0.0);
        nodes_coor[0].push_back(0.0);

        nodes_coor[1].push_back(1.0);
        nodes_coor[1].push_back(0.0);
        nodes_coor[1].push_back(0.0);

        nodes_coor[2].push_back(0.0);
        nodes_coor[2].push_back(1.0);
        nodes_coor[2].push_back(0.0);

        nodes_coor[3].push_back(1.0);
        nodes_coor[3].push_back(1.0);
        nodes_coor[3].push_back(0.0);

        nodes_coor[4].push_back(0.5);
        nodes_coor[4].push_back(0.5);
        nodes_coor[4].push_back(1.0);


        return nodes_coor;
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

    // TET_Builder no longer requires this function.
    // bool get_submesh();

};  // end class TET_test_1

}   // end namespace rtt_mc_test

#endif  // __mc_test_TET_test_1_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/TET_test_1.hh
//---------------------------------------------------------------------------//
