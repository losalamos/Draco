//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   amr_mesh/test/TestAMRMesh.hh
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the RTT_Format class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __amr_test_TestAMRMesh_hh__
#define __amr_test_TestAMRMesh_hh__

#include "UnitTestFrame/TestApp.hh"
#include "../Interface.hh"
#include "../Builder.hh"
#include "../Mesh.hh"

namespace rtt_amr_test
{
using std::map;
//===========================================================================//
/*!
 * \class TestAMRMesh
 *
 * \brief Tests the RTT_Format mesh reader class.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestAMRMesh : public rtt_UnitTestFrame::TestApp  
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    TestAMRMesh(int argc, char *argv[], std::ostream & os_in);
    // Defaulted TestAMRMesh(const TestAMRMesh &rhs);
    // Defaulted ~TestAMRMesh();

    // MANIPULATORS
    
    //Defaulted TestAMRMesh & operator=(const TestAMRMesh & rhs);

    // ACCESSORS

    std::string name() const { return "TestAMRMesh"; }

    std::string version() const;

  protected:

    std::string runTest();

  private:

    // DATA
    enum Meshes {STR1};

    // ACCESSORS
    bool check_interface(const rtt_amr::CAR_CU_Interface & interface, 
			 const rtt_meshReaders::RTT_Format & mesh, 
			 const Meshes & meshtype);


    // IMPLEMENTATION
};

} // end namespace rtt_amr_test

#endif                // __amr_test_TestAMRMesh_hh__

//---------------------------------------------------------------------------//
//                    end of amr_mesh/test/TestAMRMesh.hh
//---------------------------------------------------------------------------//
