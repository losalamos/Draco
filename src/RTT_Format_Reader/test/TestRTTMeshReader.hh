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

#include "../RTT_Mesh_Reader.hh"

namespace rtt_RTT_Mesh_Reader_test
{

typedef rtt_RTT_Format_Reader::RTT_Mesh_Reader RTT_Mesh_Reader;

enum Meshes {DEFINED};

bool check_virtual(const RTT_Mesh_Reader & mesh, const Meshes & meshtype);

} // end namespace rtt_RTT_Mesh_Reader_test

#endif                // _RTT_Format_Reader_test_TestRTT_Mesh_Reader_hh__

//---------------------------------------------------------------------------//
//             end of RTT_Format_Reader/test/TestRTTMeshReader.hh
//---------------------------------------------------------------------------//
