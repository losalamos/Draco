//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestRTT_Format.hh
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the RTT_Format class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_test_TestRTT_Format_hh__
#define __meshReaders_test_TestRTT_Format_hh__

#include "UnitTestFrame/TestApp.hh"
#include "../RTT_Format.hh"

namespace rtt_meshReaders_test
{
using std::map;
//===========================================================================//
/*!
 * \class TestRTT_Format
 *
 * \brief Tests the RTT_Format mesh reader class.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestRTT_Format : public rtt_UnitTestFrame::TestApp  
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    TestRTT_Format(int argc, char *argv[], std::ostream & os_in);
    // Defaulted TestRTT_Format(const TestRTT_Format &rhs);
    // Defaulted ~TestRTT_Format();

    // MANIPULATORS
    
    //Defaulted TestRTT_Format & operator=(const TestRTT_Format & rhs);

    // ACCESSORS

    std::string name() const { return "TestRTT_Format"; }

    std::string version() const;

  protected:

    std::string runTest();

  private:

    // DATA
    enum Meshes {DEFINED, SORTED};
    // All nested class accessor function tests with the exception of 
    // check_header and check_dims require that the Dims data has been
    // properly processed, and use the verify_Dims member function to 
    // query this data member to determine if this is true.
    map<Meshes, bool> Dims_validated;

    // ACCESSORS
    bool verify_Dims(const rtt_meshReaders::RTT_Format & mesh, 
		     const Meshes & meshtype)
    {
        // Verify that the Dims data was previously validated.
        if (!Dims_validated.count(meshtype))
	    check_dims(mesh, meshtype);
	// Return the integrity state of the Dims data.        
	return Dims_validated.find(meshtype)->second;
    }
    bool check_header(const rtt_meshReaders::RTT_Format & mesh, 
		      const Meshes & meshtype);
    bool check_dims(const rtt_meshReaders::RTT_Format & mesh, 
		    const Meshes & meshtype);
    bool check_node_flags(const rtt_meshReaders::RTT_Format & mesh, 
			  const Meshes & meshtype);
    bool check_side_flags(const rtt_meshReaders::RTT_Format & mesh, 
			  const Meshes & meshtype);
    bool check_cell_flags(const rtt_meshReaders::RTT_Format & mesh, 
			  const Meshes & meshtype);
    bool check_node_data_ids(const rtt_meshReaders::RTT_Format & mesh, 
			     const Meshes & meshtype);
    bool check_side_data_ids(const rtt_meshReaders::RTT_Format & mesh, 
			     const Meshes & meshtype);
    bool check_cell_data_ids(const rtt_meshReaders::RTT_Format & mesh, 
			     const Meshes & meshtype);
    bool check_cell_defs(const rtt_meshReaders::RTT_Format & mesh, 
			 const Meshes & meshtype);
    bool check_renumber(const rtt_meshReaders::RTT_Format & mesh, 
			const Meshes & meshtype);
    bool check_nodes(const rtt_meshReaders::RTT_Format & mesh, 
		     const Meshes & meshtype);
    bool check_sides(const rtt_meshReaders::RTT_Format & mesh, 
		     const Meshes & meshtype);
    bool check_cells(const rtt_meshReaders::RTT_Format & mesh, 
		     const Meshes & meshtype);
    bool check_node_data(const rtt_meshReaders::RTT_Format & mesh, 
			 const Meshes & meshtype);
    bool check_side_data(const rtt_meshReaders::RTT_Format & mesh, 
			 const Meshes & meshtype);
    bool check_cell_data(const rtt_meshReaders::RTT_Format & mesh, 
			 const Meshes & meshtype);
    bool check_virtual(const rtt_meshReaders::RTT_Format & mesh, 
		       const Meshes & meshtype);
    bool check_connectivity(const rtt_meshReaders::RTT_Format & mesh, 
			    const Meshes & meshtype);

    // IMPLEMENTATION
};

} // end namespace rtt_meshReaders_test

#endif                // __meshReaders_test_TestRTT_Format_hh__

//---------------------------------------------------------------------------//
//                    end of meshReaders/test/TestRTT_Format.hh
//---------------------------------------------------------------------------//
