//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTMeshReader.hh
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the RTT_Mesh_Reader class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_test_TestRTT_Mesh_Reader_hh__
#define __RTT_Format_Reader_test_TestRTT_Mesh_Reader_hh__

#include "UnitTestFrame/TestApp.hh"
#include "../RTT_Mesh_Reader.hh"

namespace rtt_RTT_Mesh_Reader_test
{
//===========================================================================//
/*!
 * \class TestRTT_Mesh_Reader
 *
 * \brief Tests the RTT_Mesh_Reader mesh reader class.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestRTT_Mesh_Reader : public rtt_UnitTestFrame::TestApp  
{

    // NESTED CLASSES AND TYPEDEFS
    typedef rtt_RTT_Format_Reader::RTT_Mesh_Reader RTT_Mesh_Reader;

    // DATA
    
  public:

    // CREATORS
    
    TestRTT_Mesh_Reader(int argc, char *argv[], std::ostream & os_in);
    // Defaulted TestRTT_Mesh_Reader(const TestRTT_Mesh_Reader &rhs);
    // Defaulted ~TestRTT_Mesh_Reader();

    // MANIPULATORS
    
    //Defaulted TestRTT_Mesh_Reader & operator=(const TestRTT_Mesh_Reader 
    //                                          & rhs);

    // ACCESSORS

    std::string name() const { return "TestRTT_Mesh_Reader"; }

    std::string version() const;

  protected:

    std::string runTest();

  private:

    // DATA
    enum Meshes {DEFINED};
    std::map<Meshes, bool> Dims_validated;

    // ACCESSORS
    bool check_virtual(const RTT_Mesh_Reader & mesh, const Meshes & meshtype);

    // IMPLEMENTATION
};

} // end namespace rtt_RTT_Mesh_Reader_test

#endif                // _RTT_Format_Reader_test_TestRTT_Mesh_Reader_hh__

//---------------------------------------------------------------------------//
//             end of RTT_Format_Reader/test/TestRTTMeshReader.hh
//---------------------------------------------------------------------------//
