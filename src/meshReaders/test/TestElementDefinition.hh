//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestElementDefinition.hh
 * \author John McGhee
 * \date   Fri Mar  3 08:41:46 2000
 * \brief  Header file for the Element_Definition class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_test_TestElementDefinition_hh__
#define __meshReaders_test_TestElementDefinition_hh__

#include "../Element_Definition.hh"

namespace rtt_meshReaders_test
{   

bool test_node(const rtt_meshReaders::Element_Definition elem_def);
bool test_bar_2(const rtt_meshReaders::Element_Definition elem_def);
bool test_bar_3(const rtt_meshReaders::Element_Definition elem_def);
bool test_tri_3(const rtt_meshReaders::Element_Definition elem_def);
bool test_tri_6(const rtt_meshReaders::Element_Definition elem_def);
bool test_quad_4(const rtt_meshReaders::Element_Definition elem_def);
bool test_quad_8(const rtt_meshReaders::Element_Definition elem_def);
bool test_quad_9(const rtt_meshReaders::Element_Definition elem_def);
bool test_tetra_4(const rtt_meshReaders::Element_Definition elem_def);
bool test_tetra_10(const rtt_meshReaders::Element_Definition elem_def);
bool test_pyra_5(const rtt_meshReaders::Element_Definition elem_def);
bool test_pyra_14(const rtt_meshReaders::Element_Definition elem_def);
bool test_penta_6(const rtt_meshReaders::Element_Definition elem_def);
bool test_penta_15(const rtt_meshReaders::Element_Definition elem_def);
bool test_penta_18(const rtt_meshReaders::Element_Definition elem_def);
bool test_hexa_8(const rtt_meshReaders::Element_Definition elem_def);
bool test_hexa_20(const rtt_meshReaders::Element_Definition elem_def);
bool test_hexa_27(const rtt_meshReaders::Element_Definition elem_def);

} // end namespace rtt_meshReaders_test

#endif           // __meshReaders_test_TestElementDefinition_hh__

//---------------------------------------------------------------------------//
//               end of meshReaders/test/TestElementDefinition.hh
//---------------------------------------------------------------------------//
