//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestHexMeshReader.hh
 * \author John McGhee
 * \date   Thu Mar  9 08:54:59 2000
 * \brief  Header file for the Hex_Mesh_Reader class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_test_TestHexMeshReader_hh__
#define __meshReaders_test_TestHexMeshReader_hh__

#include "UnitTestFrame/TestApp.hh"
#include <string>
#include <map>
#include <set>

namespace rtt_meshReaders 
{

// Forward References
class Hex_Mesh_Reader;

} //end rtt_meshReaders namespace

namespace rtt_meshReaders_test
{
 
//===========================================================================//
/*!
 * \class TestHexMeshReader
 *
 * \brief Tests the Hex_Mesh_Reader mesh reader class.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestHexMeshReader : public rtt_UnitTestFrame::TestApp  
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    TestHexMeshReader(int argc, char *argv[], std::ostream &os_in);
    // Defaulted TestHexMeshReader(const TestHexMeshReader &rhs);
    // Defaulted ~TestHexMeshReader();

    // MANIPULATORS
    
    //Defaulted TestHexMeshReader& operator=(const TestHexMeshReader &rhs);

    // ACCESSORS

    std::string name() const { return "TestHexMeshReader"; }

    std:: string version() const;

  protected:

    std::string runTest();
    
  private:

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
    bool compare_double(const double &lhs, const double &rhs);
    bool check_map(const std::map<std::string, std::set<int> >
		   &elmsets, const std::string &name, const int &begin, 
		   const int &end);

    // IMPLEMENTATION
};

} // end namespace rtt_meshReaders_test

#endif                // __meshReaders_test_TestHexMeshReader_hh__

//---------------------------------------------------------------------------//
//                    end of meshReaders/test/TestHexMeshReader.hh
//---------------------------------------------------------------------------//
