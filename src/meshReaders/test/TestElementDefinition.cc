//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestElementDefinition.cc
 * \author John McGhee
 * \date   Fri Mar  3 08:41:46 2000
 * \brief  Implements the Element_Definition unit tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestElementDefinition.hh"
#include "../Element_Definition.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include <sstream>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;


namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_meshReaders_test::TestElementDefinition;
    
    return SP<TestApp>(new TestElementDefinition(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_meshReaders_test
{

using std::string;
using rtt_meshReaders::Element_Definition;

TestElementDefinition::TestElementDefinition(int argc, char *argv[],
					       std::ostream& os_in)

    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestElementDefinition" << endl;
}

string TestElementDefinition::version() const
{
    return rtt_meshReaders::release();
}

/*!
 * \brief Tests the element definitions.
 *
 * This class builds all the elements currently defined. Currently, the
 * only real check is that an assertion doesn't fire. The output of
 * the test has been verified by hand, but a comparison with "correct"
 * results needs to be automated. (jmm, 4 Mar 00)
 *
 */
string TestElementDefinition::runTest()
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
    std::vector<int> tmp;
    cout << endl << "Building Elements for Test ---" << endl << endl;
    for (int i=0; i< type_list.size(); i++)
    {
	elem_defs.push_back( Element_Definition(type_list[i]) );
	cout << "Element Type   : " << elem_defs[i].get_type() << endl;
	cout << "Element Name   : " << elem_defs[i].get_name() << endl;
	cout << "Number of Nodes: " << elem_defs[i].get_number_of_nodes() <<
	    endl;
	cout << "Dimension      : " << elem_defs[i].get_dimension() << endl;
	cout << "Number of Sides: " << elem_defs[i].get_number_of_sides() << endl;
	cout << "Node Locations : ";
	for (int j=0; j<elem_defs[i].get_number_of_nodes(); j++)
	    cout << elem_defs[i].get_node_location(j) << " ";
	cout << endl;

	cout << "Side Types     : ";
	for (int j=0; j<elem_defs[i].get_number_of_sides(); j++)
	    cout << elem_defs[i].get_side_type(j).get_name() << " ";
	cout << endl;

	cout << "Side Nodes     : " << endl;
	for (int j=0; j<elem_defs[i].get_number_of_sides(); j++)
	{
	    tmp = elem_defs[i].get_side_nodes(j);
	    cout << "  " << "side# " << j << " -    ";
	    for (int k=0; k<tmp.size(); k++)
		cout << tmp[k] << " ";
	    cout << endl;
	}
	cout << endl;
    }
    cout << endl;
	
    pass() << "All tests passed.";
    return "All tests passed.";

//     fail() << "This test sucks!";

//     if (passed())
//     {
// 	pass() << "All tests passed.";
// 	return "All tests passed.";
//     }
//     return "Some tests failed.";
}

} // end namespace rtt_meshReaders_test


//---------------------------------------------------------------------------//
//                              end of TestElementDefinition.cc
//---------------------------------------------------------------------------//
