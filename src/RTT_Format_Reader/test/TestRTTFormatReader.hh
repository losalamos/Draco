//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTFormatReader.hh
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the RTT_Format_Reader class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __test_TestRTT_Format_Reader_hh__
#define __test_TestRTT_Format_Reader_hh__

#include "UnitTestFrame/TestApp.hh"
#include "../RTT_Format_Reader.hh"
#include <map>

namespace rtt_RTT_Format_Reader_test
{
//===========================================================================//
/*!
 * \class TestRTT_Format_Reader
 *
 * \brief Tests the RTT_Format_Reader mesh reader class.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestRTT_Format_Reader : public rtt_UnitTestFrame::TestApp  
{

    // NESTED CLASSES AND TYPEDEFS
    typedef rtt_RTT_Format_Reader::RTT_Format_Reader RTT_Format_Reader;

    // DATA
    
  public:

    // CREATORS
    
    TestRTT_Format_Reader(int argc, char *argv[], std::ostream & os_in);
    // Defaulted TestRTT_Format_Reader(const TestRTT_Format_Reader &rhs);
    // Defaulted ~TestRTT_Format_Reader();

    // MANIPULATORS
    
    //Defaulted TestRTT_Format_Reader & operator=(const TestRTT_Format_Reader 
    //                                            & rhs);

    // ACCESSORS

    std::string name() const { return "TestRTT_Format_Reader"; }

    std::string version() const;

  protected:

    std::string runTest();

  private:

    // DATA
    enum Meshes {DEFINED, SORTED, AMR};
    // All member data class accessor function tests with the exception of 
    // check_header and check_dims require that the Dims data has been
    // properly processed, and use the verify_Dims member function to 
    // query this data member to determine if this is true.
    std::map<Meshes, bool> Dims_validated;

    // ACCESSORS
    bool verify_Dims(const RTT_Format_Reader & mesh, 
		     const Meshes & meshtype)
    {
        // Verify that the Dims data was previously validated.
        if (!Dims_validated.count(meshtype))
	    check_dims(mesh, meshtype);
	// Return the integrity state of the Dims data.        
	return Dims_validated.find(meshtype)->second;
    }
    bool check_header(const RTT_Format_Reader & mesh, const Meshes & meshtype);
    bool check_dims(const RTT_Format_Reader & mesh, const Meshes & meshtype);
    bool check_node_flags(const RTT_Format_Reader & mesh, 
			  const Meshes & meshtype);
    bool check_side_flags(const RTT_Format_Reader & mesh, 
			  const Meshes & meshtype);
    bool check_cell_flags(const RTT_Format_Reader & mesh, 
			  const Meshes & meshtype);
    bool check_node_data_ids(const RTT_Format_Reader & mesh, 
			     const Meshes & meshtype);
    bool check_side_data_ids(const RTT_Format_Reader & mesh, 
			     const Meshes & meshtype);
    bool check_cell_data_ids(const RTT_Format_Reader & mesh, 
			     const Meshes & meshtype);
    bool check_cell_defs(const RTT_Format_Reader & mesh, 
			 const Meshes & meshtype);
    bool check_renumber(const RTT_Format_Reader & mesh, 
			const Meshes & meshtype);
    bool check_nodes(const RTT_Format_Reader & mesh, const Meshes & meshtype);
    bool check_sides(const RTT_Format_Reader & mesh, const Meshes & meshtype);
    bool check_cells(const RTT_Format_Reader & mesh, const Meshes & meshtype);
    bool check_node_data(const RTT_Format_Reader & mesh, 
			 const Meshes & meshtype);
    bool check_side_data(const RTT_Format_Reader & mesh, 
			 const Meshes & meshtype);
    bool check_cell_data(const RTT_Format_Reader & mesh, 
			 const Meshes & meshtype);

    // IMPLEMENTATION
};

} // end namespace rtt_RTT_Format_Reader_test

#endif                // _test_TestRTT_Format_Reader_hh__

//---------------------------------------------------------------------------//
//                    end of meshReaders/test/TestRTTFormatReader.hh
//---------------------------------------------------------------------------//
