//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstFields.cc
 * \author Thomas M. Evans
 * \date   Thu Dec 20 16:30:46 2001
 * \brief  Fields test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../OS_Mesh.hh"
#include "../OS_Builder.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc_test::Parser;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// cell-centered scalar fields - basic test of functionality

template<class MT>
void test_CCSF(SP<MT> mesh)
{
    // now build a CCSF
    typename MT::CCSF<double> field(mesh);
    if (field.empty())                     ITFAILS;
    if (field.size() != mesh->num_cells()) ITFAILS;
    if (field.get_Mesh() != *mesh)         ITFAILS;

    // lets go through this guy and fill up its stuff
    typename MT::CCSF<double>::iterator iter;
    double value = 10;
    for (iter = field.begin(); iter != field.end(); iter++)
	*iter = value++;

    // now lets check it through subscripting
    value = 10;
    for (int i = 1; i <= field.size(); i++)
	if (field(i) != value++) ITFAILS;

    // lets do some reassignment
    for (int i = 1; i <= field.size(); i++)
	field(i) *= 2;

    // lets check the reassignment by iterating through it
    value = 10;
    for (iter = field.begin(); iter != field.end(); iter++)
	if (*iter != 2 * value++) ITFAILS;

    // now we check vector assignment and iteration
    vector<double> ref(mesh->num_cells());
    for (int i = 0; i < ref.size(); i++)
	ref[i] = static_cast<double>(i) * 2;
    typename MT::CCSF<double> ccf(mesh, ref);

    value = 0;
    int index = 0;
    for (iter = ccf.begin(); iter != ccf.end(); iter++)
    {
	// check both subscripting options [] and (), notice that [] is
	// traditional C-style (0:n-1) whereas () is "cell-style" (1:N)
	if (*iter != ccf[index++]) ITFAILS;
	if (*iter != ccf(index))   ITFAILS;

	// check the actual value of the data
	if (*iter != 2 * value++)  ITFAILS;
    }
}

//---------------------------------------------------------------------------//
// cell-centered vector fields - basic test of functionality

template<class MT>
void test_CCVF(SP<MT> mesh)
{
    typename MT::CCVF<double> field(mesh);
    if (field.empty())             ITFAILS;
    if (field.size() != 2)         ITFAILS;
    if (field.get_Mesh() != *mesh) ITFAILS;

    // fill up the field
    typename MT::CCSF<double>::iterator itor;
    double value = 10;
    for (int i = 1; i <= field.size(); i++)
    {
	if (field.size(i) != mesh->num_cells()) ITFAILS;

	for (itor = field.begin(i); itor != field.end(i); itor++)
	    *itor = value * i;
    }

    // check the field
    for (int i = 1; i <= field.size(); i++)
	for (int j = 1; j <= field.size(i); j++)
	    if (field(i,j) != value * i) ITFAILS;

    // check an algorithm
    field(1, 3) = 12;
    typename MT::CCSF<double>::iterator find_itor;
    
    find_itor = find(field.begin(1), field.end(1), 12);

    if (*find_itor != 12)                ITFAILS;
    if (find_itor != field.begin(1) + 2) ITFAILS;
}

//---------------------------------------------------------------------------//
// cell-centered scalar fields - test of STL algortihm compatibility

template<class MT>
void test_CCSF_STL(SP<MT> mesh)
{
    // build a vector of info and a CCSF
    vector<double> ref(mesh->num_cells(), 10.0);
    typename MT::CCSF<double> field(mesh, ref);
    if (ref.size() != field.size()) ITFAILS;

    field(1) = 5.0;
    field(4) = 5.0;

    // make first and last iterators
    typename MT::CCSF<double>::iterator first, last, var;
    first = field.begin();
    last  = field.end();
    if (last - first != mesh->num_cells()) ITFAILS;

    // let's do some counting
    if (count(first, last, 5.0) != 2)  ITFAILS;
    if (count(first, last, 10.0) != 4) ITFAILS;

    // let's do some sorting, uses a random access iterator
    field(6) = 3.0;
    sort(first, last);
    
    if (*(first + 0) != 3.0) ITFAILS;
    if (*(first + 1) != 5.0) ITFAILS;
    if (*(first + 2) != 5.0) ITFAILS;
    for (var = first + 3; var != last; var++)
	if (*var != 10.0) ITFAILS;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is a scalar test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// 2D Mesh tests

	// build a mesh
	SP<Parser> interface(new Parser());
	OS_Builder builder(interface);
	SP<OS_Mesh> mesh = builder.build_Mesh();

	// run the tests
	test_CCSF(mesh);
	test_CCVF(mesh);
	test_CCSF_STL(mesh);
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstFields, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	C4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstFields Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstFields on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstFields.cc
//---------------------------------------------------------------------------//
