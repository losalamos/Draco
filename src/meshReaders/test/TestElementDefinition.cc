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

string TestElementDefinition::runTest()
{
    using rtt_meshReaders::Element_Definition;
    std::vector<Element_Definition> elem_defs;
    cout << endl << "Building Elements --" << endl << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::NODE));
    cout << " Built a NODE" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::BAR_2));
    cout << " Built a BAR_2" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::BAR_3));
    cout << " Built a BAR_3" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::TRI_3));
    cout << " Built a TRI_3" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::TRI_6));
    cout << " Built a TRI_6" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::QUAD_4));
    cout << " Built a QUAD_4" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::QUAD_8));
    cout << " Built a QUAD_8" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::QUAD_9));
    cout << " Built a QUAD_9" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::TETRA_4));
    cout << " Built a TETRA_4" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::TETRA_10));
    cout << " Built a TETRA_10" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::PYRA_5));
    cout << " Built a PYRA_5" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::PYRA_14));
    cout << " Built a PYRA_14" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::PENTA_6));
    cout << " Built a PENTA_6" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::PENTA_15));
    cout << " Built a PENTA_15" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::PENTA_18));
    cout << " Built a PENTA_18" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::HEXA_8));
    cout << " Built a HEXA_8" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::HEXA_20));
    cout << " Built a HEXA_20" << endl;
    elem_defs.push_back(Element_Definition(Element_Definition::HEXA_27));
    cout << " Built a HEXA_27" << endl;

    std::vector<int> tmp;
    cout << endl << "Detailed Element Descriptions --" << endl << endl;
    for (int i=0; i< elem_defs.size(); i++)
    {
	cout << "Element Name   : " << elem_defs[i].get_name() << endl;
	cout << "Number of Nodes: " << elem_defs[i].get_number_of_nodes() <<
	    endl;
	cout << "Dimension      : " << elem_defs[i].get_dimension() << endl;
	cout << "Number of Sides: " << elem_defs[i].get_nsides() << endl;
	cout << "Node Locations : ";
	for (int j=0; j<elem_defs[i].get_number_of_nodes(); j++)
	    cout << elem_defs[i].get_node_loc(j) << " ";
	cout << endl;

	cout << "Side Types     : ";
	for (int j=0; j<elem_defs[i].get_nsides(); j++)
	    cout << elem_defs[i].get_side_type(j).get_name() << " ";
	cout << endl;

	cout << "Side Nodes     : " << endl;
	for (int j=0; j<elem_defs[i].get_nsides(); j++)
	{
	    tmp = elem_defs[i].get_side_nodes(j);
	    cout << "  ";
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
