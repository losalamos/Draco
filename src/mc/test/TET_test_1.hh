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

namespace rtt_mc_test
{

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
//  0) original : Committed 2000_02_09
//___________________________________________________________________________//

class TET_test_1
{

 public:

    // Default constructor.  No need to initialize anything.
    TET_test_1() { }

    std::string get_coord_system() { return std::string("XYZ"); }

    std::vector< std::vector<double> > get_nodes()
    {
        std::vector< std::vector<double> > nodes;
        nodes.resize(3);

        nodes[0].push_back(0.0);
        nodes[0].push_back(1.0);
        nodes[0].push_back(0.0);
        nodes[0].push_back(1.0);
        nodes[0].push_back(0.5);

        nodes[1].push_back(0.0);
        nodes[1].push_back(0.0);
        nodes[1].push_back(1.0);
        nodes[1].push_back(1.0);
        nodes[1].push_back(0.5);

        nodes[2].push_back(0.0);
        nodes[2].push_back(0.0);
        nodes[2].push_back(0.0);
        nodes[2].push_back(0.0);
        nodes[2].push_back(1.0);

        return nodes;
    }

    std::vector<int> get_parent()
    {
        std::vector<int> parent;

        parent.push_back(1);
        parent.push_back(2);
        parent.push_back(3);
        parent.push_back(4);
        parent.push_back(5);

        return parent;
    }

    std::vector< std::vector<int> > get_cells_nodes()
    {
        std::vector< std::vector<int> > cells_nodes;
        cells_nodes.resize(2);

        cells_nodes[0].push_back(1);
        cells_nodes[0].push_back(2);
        cells_nodes[0].push_back(3);
        cells_nodes[0].push_back(5);

        cells_nodes[1].push_back(3);
        cells_nodes[1].push_back(2);
        cells_nodes[1].push_back(5);
        cells_nodes[1].push_back(4);

        return cells_nodes;
    }

    bool get_submesh() { return false; }

};  // end class TET_test_1

}   // end namespace rtt_mc_test

#endif  // __mc_test_TET_test_1_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/TET_test_1.hh
//---------------------------------------------------------------------------//
