//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/test/tstEnsightTranslator.cc
 * \author Thomas M. Evans
 * \date   Mon Jan 24 11:12:59 2000
 * \brief  Ensight_Translator test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "viz_test.hh"
#include "../Ensight_Translator.hh"
#include "../Release.hh"
#include "ds++/Assert.hh"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h> 

using namespace std;
using rtt_viz::Ensight_Translator;

//---------------------------------------------------------------------------//

void ensight_dump_test(const bool binary)
{
    // dimensions
    int ncells   = 27; 
    int nvert    = 64; 
    int ndim     = 3;
    int ndata    = 2;
    int nhexvert = 8;
    int nrgn     = 2;
    
    // do an Ensight Dump
    vector<vector<int> >    ipar(ncells, vector<int>(nhexvert));
    vector<vector<double> > ens_vrtx_data(nvert, vector<double>(ndata, 5.0)); 
    vector<vector<double> > ens_cell_data(ncells, vector<double>(ndata, 10.));   
    vector<vector<double> > pt_coor(nvert, vector<double>(ndim));
  
    vector<int>    iel_type(ncells, 
    			    rtt_viz::eight_node_hexahedron);
    vector<int>    rgn_index(ncells, 1);
    vector<string> ens_vdata_names(ndata, "Temperatures");
    vector<string> ens_cdata_names(ndata, "Velocity");
    vector<string> rgn_name(nrgn, "RGN_A");
    vector<int>    rgn_data(nrgn, 1);

    // set region stuff
    rgn_name[1] = "RGN_B";
    rgn_data[1] = 2;
    for (int i = 1; i < 5; i++)
	rgn_index[i] = 2;
    rgn_index[14] = 2;
    rgn_index[15] = 2;
    rgn_index[21] = 2;
    ens_vdata_names[1] = "Densities";
    ens_cdata_names[1] = "Pressure";

    string prefix   = "testproblem";
    if ( binary )
	prefix += "_binary";
    
    int icycle      = 1;
    double time     = .01;
    double dt       = .01;
    string gd_wpath = ".";

    // make data
    for (int i = 0; i < ndata; i++)
    {
	// cell data
	for (int cell = 0; cell < ncells; cell++)
	    ens_cell_data[cell][i] = 1 + cell;

	// vrtx data
	for (int v = 0; v < nvert; v++)
	    ens_vrtx_data[v][i] = 1 + v;
    }

    // read cell data
    ifstream input("cell_data");
    for (int i = 0; i < pt_coor.size(); i++)
	for (int j = 0; j < pt_coor[i].size(); j++)
	    input >> pt_coor[i][j];
    for (int i = 0; i < ipar.size(); i++)
	for (int j = 0; j < ipar[i].size(); j++)
	    input >> ipar[i][j];

    const bool static_geom = false;

    // build an Ensight_Translator (make sure it overwrites any existing
    // stuff) 
    Ensight_Translator translator(prefix, gd_wpath, ens_vdata_names,
				  ens_cdata_names, true, static_geom,
				  binary); 

    translator.ensight_dump(icycle, time, dt,
			    ipar, iel_type, rgn_index, pt_coor,
			    ens_vrtx_data, ens_cell_data,
			    rgn_data, rgn_name);

    vector<double> dump_times = translator.get_dump_times();
    if (dump_times.size() != 1) ITFAILS;
    if (dump_times[0] != .01)   ITFAILS;

    // build another ensight translator; this should overwrite the existing
    // directories
    Ensight_Translator translator2(prefix, gd_wpath, ens_vdata_names,
				   ens_cdata_names, true, static_geom,
				   binary); 
    
    translator2.ensight_dump(icycle, time, dt,
			     ipar, iel_type, rgn_index, pt_coor,
			     ens_vrtx_data, ens_cell_data,
			     rgn_data, rgn_name);

    // build another ensight translator from the existing dump times list;
    // thus we will not overwrite the existing directories

    Ensight_Translator translator3(prefix, gd_wpath, ens_vdata_names,
 				   ens_cdata_names, false, static_geom,
				   binary); 
    
    // now add another dump to the existing data
    translator3.ensight_dump(2, .05, dt,
			     ipar, iel_type, rgn_index, pt_coor,
			     ens_vrtx_data, ens_cell_data,
			     rgn_data, rgn_name);    

    // make yet a fourth translator that will append
    Ensight_Translator translator4(prefix, gd_wpath, ens_vdata_names,
				   ens_cdata_names, false, static_geom,
				   binary); 
    
    // add yet another dump to the existing data
    translator4.ensight_dump(3, .10, dt,
			     ipar, iel_type, rgn_index, pt_coor,
			     ens_vrtx_data, ens_cell_data,
			     rgn_data, rgn_name);    
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_viz::release() << endl; 
	    return 0;
	}

    try
    {
	// tests
	ensight_dump_test(false); // ascii dump
	
	// run python diff scrips (only works for ascii)
	system("python ./tstEnsight_Diff.py");

	ensight_dump_test(true); // binary dump

	// ... there's no check for binary, yet.
    }
    catch(rtt_dsxx::assertion &ass)
    {
	cout << "Dumbass you screwed up on " << ass.what() << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "**********************************************" << endl;
    if (rtt_viz_test::passed) 
    {
        cout << "**** Ensight_Translator Self Test: PASSED ****" << endl;
    }
    cout <<     "**********************************************" << endl;
    cout << endl;
    
    cout << "Done testing Ensight_Translator." << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstEnsight_Translator.cc
//---------------------------------------------------------------------------//
