//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestElementDefinition.cc
 * \author Thomas M. Evans
 * \date   Tue Mar 26 16:06:55 2002
 * \brief  Test Element Definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "meshReaders_test.hh"
#include "TestElementDefinition.hh"
#include "../Element_Definition.hh"
#include "../Release.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void runTest()
{
    using rtt_meshReaders::Element_Definition;

    std::vector<Element_Definition::Element_Type> type_list;
    type_list.push_back(Element_Definition::NODE);
    type_list.push_back(Element_Definition::BAR_2);
    type_list.push_back(Element_Definition::BAR_3);
    type_list.push_back(Element_Definition::TRI_3);
    type_list.push_back(Element_Definition::TRI_6);
    type_list.push_back(Element_Definition::QUAD_4);
    type_list.push_back(Element_Definition::QUAD_8);
    type_list.push_back(Element_Definition::QUAD_9);
    type_list.push_back(Element_Definition::TETRA_4);
    type_list.push_back(Element_Definition::TETRA_10);
    type_list.push_back(Element_Definition::PYRA_5);
    type_list.push_back(Element_Definition::PYRA_14);
    type_list.push_back(Element_Definition::PENTA_6);
    type_list.push_back(Element_Definition::PENTA_15);
    type_list.push_back(Element_Definition::PENTA_18);
    type_list.push_back(Element_Definition::HEXA_8);
    type_list.push_back(Element_Definition::HEXA_20);
    type_list.push_back(Element_Definition::HEXA_27);

    std::vector<Element_Definition> elem_defs;
    cout << endl << "Building Elements for Test ---" << endl << endl;
    for (int i=0; i< type_list.size(); i++)
    {
	elem_defs.push_back( Element_Definition(type_list[i]) );
	cout << elem_defs[i];
    }
    cout << endl;

    cout << "Checking Elements ---" << endl << endl;

    rtt_meshReaders_test::test_node(elem_defs[0]);
    rtt_meshReaders_test::test_bar_2(elem_defs[1]);
    rtt_meshReaders_test::test_bar_3(elem_defs[2]);
    rtt_meshReaders_test::test_tri_3(elem_defs[3]);
    rtt_meshReaders_test::test_tri_6(elem_defs[4]);
    rtt_meshReaders_test::test_quad_4(elem_defs[5]);
    rtt_meshReaders_test::test_quad_8(elem_defs[6]);
    rtt_meshReaders_test::test_quad_9(elem_defs[7]);
    rtt_meshReaders_test::test_tetra_4(elem_defs[8]);
    rtt_meshReaders_test::test_tetra_10(elem_defs[9]);
    rtt_meshReaders_test::test_pyra_5(elem_defs[10]);
    rtt_meshReaders_test::test_pyra_14(elem_defs[11]);
    rtt_meshReaders_test::test_penta_6(elem_defs[12]);
    rtt_meshReaders_test::test_penta_15(elem_defs[13]);
    rtt_meshReaders_test::test_penta_18(elem_defs[14]);
    rtt_meshReaders_test::test_hexa_8(elem_defs[15]);
    rtt_meshReaders_test::test_hexa_20(elem_defs[16]);
    rtt_meshReaders_test::test_hexa_27(elem_defs[17]);

    if (rtt_meshReaders_test::passed)
    {
	PASSMSG("All tests passed.");
    }
    else
    {
	FAILMSG("Some tests failed.");
    }
}

//---------------------------------------------------------------------------//

namespace rtt_meshReaders_test
{

bool test_node(const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the NODE element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="NODE";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::NODE;
    ldum = ldum && elem_def.get_number_of_nodes() == 1;
    ldum = ldum && elem_def.get_dimension() == 0;
    ldum = ldum && elem_def.get_node_location(0) ==
	Element_Definition::CORNER;
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_bar_2(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the BAR_2 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="BAR_2";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::BAR_2 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 2;
    ldum = ldum && elem_def.get_dimension() == 1;
    ldum = ldum && elem_def.get_number_of_sides() == 2;
    for (int j=0; j<2; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::NODE;
    }
    const int size = 1;
    int s0[size] = {0};
    int s1[size] = {1};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_bar_3(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the BAR_3 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="BAR_3";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::BAR_3 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 3;
    ldum = ldum && elem_def.get_dimension() == 1;
    ldum = ldum && elem_def.get_number_of_sides() == 2;
    for (int j=0; j<2; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::NODE;
    }
    ldum = ldum && elem_def.get_node_location(2) == 
	Element_Definition::EDGE;
    const int size = 1;
    int s0[size] = {0};
    int s1[size] = {1};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    if (ldum)
    {
	ostringstream message; 
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_tri_3(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the TRI_3 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="TRI_3";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::TRI_3 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 3;
    ldum = ldum && elem_def.get_dimension() == 2;
    ldum = ldum && elem_def.get_number_of_sides() == 3;
    for (int j=0; j<3; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::BAR_2;
    }
    const int size = 2;
    int s0[size] = {0,1};
    int s1[size] = {1,2};
    int s2[size] = {2,0};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    if (ldum)
    {
	ostringstream message; 
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_tri_6(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="TRI_6";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::TRI_6 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 6;
    ldum = ldum && elem_def.get_dimension() == 2;
    ldum = ldum && elem_def.get_number_of_sides() == 3;
    for (int j=0; j<3; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::BAR_3;
    }
    for (int j=3; j<6; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    const int size = 3;
    int s0[size] = {0,1,3};
    int s1[size] = {1,2,4};
    int s2[size] = {2,0,5};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_quad_4(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the QUAD_4 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="QUAD_4";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::QUAD_4 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 4;
    ldum = ldum && elem_def.get_dimension() == 2;
    ldum = ldum && elem_def.get_number_of_sides() == 4;
    for (int j=0; j<4; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::BAR_2;
    }
    const int size = 2;
    int s0[size] = {0,1};
    int s1[size] = {1,2};
    int s2[size] = {2,3};
    int s3[size] = {3,0};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+size);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_quad_8(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the QUAD_8 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="QUAD_8";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::QUAD_8 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 8;
    ldum = ldum && elem_def.get_dimension() == 2;
    ldum = ldum && elem_def.get_number_of_sides() == 4;
    for (int j=0; j<4; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::BAR_3;
    }
    for (int j=4; j<8; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    const int size = 3;
    int s0[size] = {0,1,4};
    int s1[size] = {1,2,5};
    int s2[size] = {2,3,6};
    int s3[size] = {3,0,7};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+size);
    if (ldum)
    {
	ostringstream message; 
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}   

bool test_quad_9(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the QUAD_9 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="QUAD_9";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::QUAD_9 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 9;
    ldum = ldum && elem_def.get_dimension() == 2;
    ldum = ldum && elem_def.get_number_of_sides() == 4;
    for (int j=0; j<4; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::BAR_3;
    }
    for (int j=4; j<8; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    ldum = ldum && elem_def.get_node_location(8) == 
	Element_Definition::FACE;
    const int size = 3;
    int s0[size] = {0,1,4};
    int s1[size] = {1,2,5};
    int s2[size] = {2,3,6};
    int s3[size] = {3,0,7};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+size);
    if (ldum)
    {
	ostringstream message; 
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_tetra_4(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the TETRA_4 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="TETRA_4";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::TETRA_4 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 4;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 4;
    for (int j=0; j<4; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::TRI_3;
    }
    const int size = 3;
    int s0[size] = {0,2,1};
    int s1[size] = {0,1,3};
    int s2[size] = {1,2,3};
    int s3[size] = {2,0,3};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+size);
    if (ldum)
    {
	ostringstream message; 
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_tetra_10(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the TETRA_10 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="TETRA_10";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::TETRA_10 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 10;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 4;
    for (int j=0; j<4; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::TRI_6;
    }
    for (int j=4; j<10; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    const int size = 6;
    int s0[size] = {0,2,1,6,5,4};
    int s1[size] = {0,1,3,4,8,7};
    int s2[size] = {1,2,3,5,9,8};
    int s3[size] = {2,0,3,6,7,9};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+size);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_pyra_5(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the PYRA_5 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="PYRA_5";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::PYRA_5 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 5;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 5;
    ldum = ldum && elem_def.get_node_location(0) == 
	Element_Definition::CORNER;
    ldum = ldum && elem_def.get_side_type(0).get_type() == 
	Element_Definition::QUAD_4;
    for (int j=1; j<5; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::TRI_3;
    }
    const int sizeq = 4;
    int s0[sizeq] = {0,3,2,1};
    const int sizet = 3;
    int s1[sizet] = {0,1,4};
    int s2[sizet] = {1,2,4};
    int s3[sizet] = {2,3,4};
    int s4[sizet] = {3,0,4};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+sizeq);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+sizet);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+sizet);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+sizet);
    ldum = ldum && elem_def.get_side_nodes(4) == 
	std::vector<int>(s4,s4+sizet);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_pyra_14(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the PYRA_14 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="PYRA_14";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::PYRA_14 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 14;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 5;
    ldum = ldum && elem_def.get_node_location(0) == 
	Element_Definition::CORNER;
    ldum = ldum && elem_def.get_side_type(0).get_type() == 
	Element_Definition::QUAD_8;
    for (int j=1; j<5; ++j)
    {
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::TRI_6;
    }
    for (int j=5; j<13; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    ldum = ldum && elem_def.get_node_location(13) == 
	Element_Definition::CELL;
    const int sizeq = 8;
    int s0[sizeq] = {0,3,2,1,8,7,6,5};
    const int sizet = 6;
    int s1[sizet] = {0,1,4,5,10,9};
    int s2[sizet] = {1,2,4,6,11,10};
    int s3[sizet] = {2,3,4,7,12,11};
    int s4[sizet] = {3,0,4,8,9,12};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+sizeq);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+sizet);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+sizet);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+sizet);
    ldum = ldum && elem_def.get_side_nodes(4) == 
	std::vector<int>(s4,s4+sizet);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_penta_6(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the PENTA_6 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="PENTA_6";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::PENTA_6 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 6;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 5;
    for (int j=0; j<6; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
    for (int j=0; j<3; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::QUAD_4;
    for (int j=3; j<5; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::TRI_3;
    const int sizeq = 4;
    int s0[sizeq] = {0,1,4,3};
    int s1[sizeq] = {1,2,5,4};
    int s2[sizeq] = {2,0,3,5};
    const int sizet = 3;
    int s3[sizet] = {0,2,1};
    int s4[sizet] = {3,4,5};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+sizeq);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+sizeq);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+sizeq);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+sizet);
    ldum = ldum && elem_def.get_side_nodes(4) == 
	std::vector<int>(s4,s4+sizet);
    if (ldum)
    {
	ostringstream message; 
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_penta_15(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the PENTA_15 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="PENTA_15";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::PENTA_15 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 15;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 5;
    for (int j=0; j<6; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
    for (int j=6; j<15; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    for (int j=0; j<3; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::QUAD_8;
    for (int j=3; j<5; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::TRI_6;
    const int sizeq = 8;
    int s0[sizeq] = {0,1,4,3,6,10,12,9};
    int s1[sizeq] = {1,2,5,4,7,11,13,10};
    int s2[sizeq] = {2,0,3,5,8,9,14,11};
    const int sizet = 6;
    int s3[sizet] = {0,2,1,8,7,6};
    int s4[sizet] = {3,4,5,12,13,14};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+sizeq);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+sizeq);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+sizeq);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+sizet);
    ldum = ldum && elem_def.get_side_nodes(4) == 
	std::vector<int>(s4,s4+sizet);
    if (ldum)
    {
	ostringstream message; 
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_penta_18(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the PENTA_18 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="PENTA_18";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::PENTA_18 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 18;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 5;
    for (int j=0; j<6; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
    for (int j=6; j<15; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    for (int j=15; j<18; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::FACE;
    for (int j=0; j<3; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::QUAD_9;
    for (int j=3; j<5; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::TRI_6;
    const int sizeq = 9;
    int s0[sizeq] = {0,1,4,3,6,10,12,9,15};
    int s1[sizeq] = {1,2,5,4,7,11,13,10,16};
    int s2[sizeq] = {2,0,3,5,8,9,14,11,17};
    const int sizet = 6;
    int s3[sizet] = {0,2,1,8,7,6};
    int s4[sizet] = {3,4,5,12,13,14};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+sizeq);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+sizeq);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+sizeq);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+sizet);
    ldum = ldum && elem_def.get_side_nodes(4) == 
	std::vector<int>(s4,s4+sizet);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_hexa_8(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the HEXA_8 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="HEXA_8";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::HEXA_8 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 8;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 6;
    for (int j=0; j<8; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
    for (int j=0; j<6; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::QUAD_4;
    const int size = 4;
    int s0[size] = {0,3,2,1};
    int s1[size] = {0,1,5,4};
    int s2[size] = {1,2,6,5};
    int s3[size] = {2,3,7,6};
    int s4[size] = {0,4,7,3};
    int s5[size] = {4,5,6,7};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+size);
    ldum = ldum && elem_def.get_side_nodes(4) == 
	std::vector<int>(s4,s4+size);
    ldum = ldum && elem_def.get_side_nodes(5) == 
	std::vector<int>(s5,s5+size);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_hexa_20(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the HEXA_20 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="HEXA_20";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::HEXA_20 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 20;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 6;
    for (int j=0; j<8; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
    for (int j=8; j<20; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    for (int j=0; j<6; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::QUAD_8;
    const int size = 8;
    int s0[size] = {0,3,2,1,11,10,9,8};
    int s1[size] = {0,1,5,4,8,13,16,12};
    int s2[size] = {1,2,6,5,9,14,17,13};
    int s3[size] = {2,3,7,6,10,15,18,14};
    int s4[size] = {0,4,7,3,12,19,15,11};
    int s5[size] = {4,5,6,7,16,17,18,19};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+size);
    ldum = ldum && elem_def.get_side_nodes(4) == 
	std::vector<int>(s4,s4+size);
    ldum = ldum && elem_def.get_side_nodes(5) == 
	std::vector<int>(s5,s5+size);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

bool test_hexa_27(
    const rtt_meshReaders::Element_Definition elem_def)
{
    // Test the HEXA_27 element.
    using rtt_meshReaders::Element_Definition;
    std::string ename="HEXA_27";
    bool ldum = elem_def.get_name() == ename;
    ldum = ldum && elem_def.get_type() == Element_Definition::HEXA_27 ;
    ldum = ldum && elem_def.get_number_of_nodes() == 27;
    ldum = ldum && elem_def.get_dimension() == 3;
    ldum = ldum && elem_def.get_number_of_sides() == 6;
    for (int j=0; j<8; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::CORNER;
    for (int j=8; j<20; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::EDGE;
    for (int j=20; j<26; ++j)
	ldum = ldum && elem_def.get_node_location(j) == 
	    Element_Definition::FACE;
    ldum = ldum && elem_def.get_node_location(26) == 
	Element_Definition::CELL;
    for (int j=0; j<6; ++j)
	ldum = ldum && elem_def.get_side_type(j).get_type() == 
	    Element_Definition::QUAD_9;
    const int size = 9;
    int s0[size] = {0,3,2,1,11,10,9,8,20};
    int s1[size] = {0,1,5,4,8,13,16,12,21};
    int s2[size] = {1,2,6,5,9,14,17,13,22};
    int s3[size] = {2,3,7,6,10,15,18,14,23};
    int s4[size] = {0,4,7,3,12,19,15,11,24};
    int s5[size] = {4,5,6,7,16,17,18,19,25};
    ldum = ldum && elem_def.get_side_nodes(0) == 
	std::vector<int>(s0,s0+size);
    ldum = ldum && elem_def.get_side_nodes(1) == 
	std::vector<int>(s1,s1+size);
    ldum = ldum && elem_def.get_side_nodes(2) == 
	std::vector<int>(s2,s2+size);
    ldum = ldum && elem_def.get_side_nodes(3) == 
	std::vector<int>(s3,s3+size);
    ldum = ldum && elem_def.get_side_nodes(4) == 
	std::vector<int>(s4,s4+size);
    ldum = ldum && elem_def.get_side_nodes(5) == 
	std::vector<int>(s5,s5+size);
    if (ldum) 
    {
	ostringstream message;
	message << ename << " Element OK." << endl;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Error in " << ename << " Element." << endl;
	FAILMSG(message.str());
    }
    return ldum;
}

}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_meshReaders::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	runTest();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing TestElementDefinition, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_meshReaders_test::passed) 
    {
        cout << "**** TestElementDefinition Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing TestElementDefinition." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of TestElementDefinition.cc
//---------------------------------------------------------------------------//
