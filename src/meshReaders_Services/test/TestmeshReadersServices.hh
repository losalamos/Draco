//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders_Services/test/TestmeshReadersServices.hh
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the meshReaders_Services class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __test_TestmeshReadersServices_hh__
#define __test_TestmeshReadersServices_hh__

#include "UnitTestFrame/TestApp.hh"
#include "meshReaders/RTT_Mesh_Reader.hh"
#include "../Connect.hh"
#include "ds++/SP.hh"
#include <map>

namespace rtt_meshReaders_Services_test
{
//===========================================================================//
/*!
 * \class TestmeshReaders_Services
 *
 * \brief Tests the meshReaders_Services mesh reader class.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestmeshReaders_Services : public rtt_UnitTestFrame::TestApp  
{
    // NESTED CLASSES AND TYPEDEFS
    typedef rtt_meshReaders::RTT_Mesh_Reader RTT_Mesh_Reader;
    typedef rtt_meshReaders_Services::Connect Connect;
    typedef std::vector<int> vector_int;
    typedef std::vector<std::vector<int> > vector_vector_int;
    typedef std::vector<std::vector<std::vector<int> > > 
        vector_vector_vector_int;
    typedef std::vector<double> vector_dbl;
    typedef std::vector<std::vector<double> > vector_vector_dbl;

    // DATA
    
  public:

    // CREATORS
    
    TestmeshReaders_Services(int argc, char *argv[], std::ostream & os_in);
    // Defaulted TestmeshReaders_Services(const TestmeshReaders_Services &rhs);
    // Defaulted ~TestmeshReaders_Services();

    // MANIPULATORS

    //Defaulted TestmeshReaders_Services & operator=(const 
    //    TestmeshReaders_Services & rhs);

    // ACCESSORS

    std::string name() const { return "TestmeshReaders_Services"; }

    std::string version() const;

  protected:

    std::string runTest();

  private:

    // DATA
    enum Meshes {DEFINED, AMR};
    std::vector<std::string> bndry_flags;
    std::vector<int> flag_numbs;
    rtt_dsxx::SP<RTT_Mesh_Reader> mesh;
    rtt_dsxx::SP<Connect> connect;

    // ACCESSORS
    bool check_connect(rtt_dsxx::SP<RTT_Mesh_Reader> mesh_, 
		       rtt_dsxx::SP<Connect> connect_, 
		       const Meshes & meshtype_);

    // IMPLEMENTATION
};

} // end namespace rtt_meshReaders_Services_test

#endif                // _test_TestmeshReadersServices_hh__

//---------------------------------------------------------------------------//
//               end of meshReaders/test/TestmeshReadersServices.hh
//---------------------------------------------------------------------------//
