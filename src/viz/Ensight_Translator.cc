//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Ensight_Translator.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:36:10 2000
 * \brief  Ensight_Translator implementation file (non-templated code).
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Ensight_Translator.hh"

namespace rtt_viz
{

using std::endl;
using std::setw;
using std::ios;
using std::ofstream;
using std::string;
using std::setiosflags;

//---------------------------------------------------------------------------//
// CREATEFILENAMES PRIVATE
//---------------------------------------------------------------------------//
/*!
 * \brief Creates some of the file prefixes and filenames for ensight dump.
 *
 * \param prefix std_string giving the name of the problem
 * \param gd_wpath directory where dumps are stored
 */
void Ensight_Translator::createFilenames(const std_string &prefix,
					 const std_string &gd_wpath)
{
    // ensight directory name
    ens_prefix = gd_wpath + "/" + prefix + "_ensight";

    // case file name
    case_filename = ens_prefix + "/" + prefix + ".case";
}

//---------------------------------------------------------------------------//
// INITIALIZE PRIVATE
//---------------------------------------------------------------------------//
/*!
 * \brief Common initializer for constructors.
 *
 * \param graphics_continue If true, use existing ensight directory.
 * If false, create or wipe out the existing directory.
 */
void Ensight_Translator::initialize(const bool graphics_continue)
{
    using std::strerror;

    num_ensight_cell_types = 15;

    // Assign values to ensight_cell_names. These are
    // the official "Ensight" names that must be used in 
    // the Ensight file.
    ensight_cell_names.resize(num_ensight_cell_types);
    ensight_cell_names[0]  = "point";
    ensight_cell_names[1]  = "bar2";
    ensight_cell_names[2]  = "bar3";
    ensight_cell_names[3]  = "tria3";
    ensight_cell_names[4]  = "tria6";
    ensight_cell_names[5]  = "quad4";
    ensight_cell_names[6]  = "quad8";
    ensight_cell_names[7]  = "tetra4";
    ensight_cell_names[8]  = "tetra10";
    ensight_cell_names[9]  = "pyramid5";
    ensight_cell_names[10] = "pyramid13";
    ensight_cell_names[11] = "hexa8";
    ensight_cell_names[12] = "hexa20";
    ensight_cell_names[13] = "penta6";
    ensight_cell_names[14] = "penta15";

    // Assign values to vrtx_count, the number of vertices 
    // in a cell.
    vrtx_cnt.resize(num_ensight_cell_types);
    vrtx_cnt[0]  = 1;
    vrtx_cnt[1]  = 2;
    vrtx_cnt[2]  = 3;
    vrtx_cnt[3]  = 3;
    vrtx_cnt[4]  = 6;
    vrtx_cnt[5]  = 4;
    vrtx_cnt[6]  = 8;
    vrtx_cnt[7]  = 4;
    vrtx_cnt[8]  = 10;
    vrtx_cnt[9]  = 5;
    vrtx_cnt[10] = 13;
    vrtx_cnt[11] = 8;
    vrtx_cnt[12] = 20;
    vrtx_cnt[13] = 6;
    vrtx_cnt[14] = 15;

    // Assign values to cell_type_index. The user will
    // use these to identify cell types.
    cell_type_index.resize(num_ensight_cell_types);
    cell_type_index[0]  = point;
    cell_type_index[1]  = two_node_bar;
    cell_type_index[2]  = three_node_bar;
    cell_type_index[3]  = three_node_triangle;
    cell_type_index[4]  = six_node_triangle;
    cell_type_index[5]  = four_node_quadrangle;
    cell_type_index[6]  = eight_node_quadrangle;
    cell_type_index[7]  = four_node_tetrahedron;
    cell_type_index[8]  = ten_node_tetrahedron;
    cell_type_index[9]  = five_node_pyramid;
    cell_type_index[10] = thirteen_node_pyramid;
    cell_type_index[11] = eight_node_hexahedron;
    cell_type_index[12] = twenty_node_hexahedron;
    cell_type_index[13] = six_node_wedge;
    cell_type_index[14] = fifteen_node_wedge;

    // build the ensight directory if this is not a continuation
    if (!graphics_continue)
    { 
	// remove old ensight directory
	std::ostringstream rm_ensight;
	rm_ensight << "rm -rf " << ens_prefix;
	system(rm_ensight.str().c_str());

	int err = mkdir(ens_prefix.c_str(), ENSIGHT_DIR_MODE);
	if (err == -1)
	{
	    std::ostringstream dir_error;
	    dir_error << "Error opening ensight directory: " 
		      << strerror(errno);
	    Insist (0,  dir_error.str().c_str());
	}
    }
       
    // Check to make sure the variable names are of acceptable length 
    // and contain no forbidden characters. Ensight prohibits the
    // characters "( ) [ ] + - @ ! # * ^ $ / space", and requires
    // the names be 19 characters or less. Moreover, since these
    // names will be used to label output and to create directories,
    // the names should also be unique.

    int nens_vdata = ens_vdata_names.size();
    int nens_cdata = ens_cdata_names.size();

    typedef std::vector<sf_string::iterator> SFS_iter_vec;
    // Create a name list for testing.
    sf_string name_tmp(nens_vdata+nens_cdata);
    for (int i=0; i < nens_vdata; i++ )
	name_tmp[i] = ens_vdata_names[i];
    for (int i=0; i < nens_cdata; i++ )
	name_tmp[i+nens_vdata] = ens_cdata_names[i];
    // Check for name lengths out of limits
    {
	int low = 1;
	int high= 19;
	SFS_iter_vec result =
	    rtt_dsxx::check_string_lengths(name_tmp.begin(), 
					   name_tmp.end(), low, high);
	if (result.size() != 0) 
	{
	    std::cerr << "*** Error in variable name(s) -" << std::endl;
	    for (int i=0; i<result.size(); i++)
		std::cerr << "Size of name is not in allowable range: \"" 
			  << *result[i] << "\"" << std::endl;
	    std::cerr << "Name lengths must be greater than " << low 
		      << " and less than " << high << "." << std::endl;
	    Insist (0, "Ensight variable name length out of limits!");
	}
    }
    // Check for bad characters.
    {
	std::string bad_chars = "()[]+-@!#*^$/ ";
	SFS_iter_vec result =
	    rtt_dsxx::check_string_chars(name_tmp.begin(), 
					 name_tmp.end(), bad_chars);
	if (result.size() != 0) 
	{
	    std::cerr << "*** Error in variable name(s) -" << std::endl;
	    for (int i=0; i<result.size(); i++) 
	        std::cerr << "Found disallowed character(s) in name: \"" 
	                  << *result[i] << "\"" << std::endl;
	    std::cerr << "The following characters are forbidden:" << 
		std::endl << " \"" << bad_chars << "\"," << 
		" as well as any white-space characters." << std::endl;
	    Insist (0, "Found illegal character in ensight variable names!");
	}
    }
    // Check for non-unique names
    {
	SFS_iter_vec result =
	    rtt_dsxx::check_strings_unique(name_tmp.begin(), name_tmp.end());
	if (result.size() != 0) 
	{
	    std::cerr << "*** Error in variable name(s) -" << std::endl;
	    for (int i=0; i<result.size(); i++)
		std::cerr << "Duplicate name found: \"" 
			  << *result[i] << "\"" << std::endl;
	    std::cerr << "All variable names must be unique!" << std::endl;
	    Insist (0, "Duplicate ensight variable names found!");
	}
    }

    // calculate and make the geometry directory if this is not a
    // continuation
    geo_dir = ens_prefix + "/geo";
    if (!graphics_continue)
	mkdir(geo_dir.c_str(), ENSIGHT_DIR_MODE);

    // make data directory names and directories
    vdata_dirs.resize(ens_vdata_names.size());
    cdata_dirs.resize(ens_cdata_names.size());
    for (int i = 0; i < ens_vdata_names.size(); i++)
    {
	vdata_dirs[i] = ens_prefix + "/" + ens_vdata_names[i];
	
	// if this is not a continuation make the directory
	if (!graphics_continue)
	    mkdir(vdata_dirs[i].c_str(), ENSIGHT_DIR_MODE);
    }
    for (int i = 0; i < ens_cdata_names.size(); i++)
    {
	cdata_dirs[i] = ens_prefix + "/" + ens_cdata_names[i];

	// if this is not a continuation make the directory
	if (!graphics_continue)
	    mkdir(cdata_dirs[i].c_str(), ENSIGHT_DIR_MODE);
    }   
}

//---------------------------------------------------------------------------//
// WRITE ENSIGHT CASE FILE
//---------------------------------------------------------------------------//
/*!
 * \brief Write out case file.
 */
void Ensight_Translator::ensight_case(const double)
{
    // create the case file name (directory already created)
    const char *filename = case_filename.c_str();
    
    // open the case file
    ofstream caseout(filename);

    // write the format header
    caseout << "FORMAT" << endl;
    caseout << "type: ensight" << endl << endl;

    // write the geometry file block
    caseout << "GEOMETRY" << endl;
    caseout << "model: 1   " << "./geo/data";
    if ( ! static_geom ) {
	caseout << ".****";
    }
    caseout << endl << endl;

    // write the variable block header
    caseout << "VARIABLE" << endl;

    // write the pointer to the node variables
    for (int i = 0; i < ens_vdata_names.size(); i++)
	caseout << " scalar per node:    1  " << setw(19)
		<< setiosflags(ios::left)
		<< ens_vdata_names[i] << setw(4) << " "
		<< "./" << ens_vdata_names[i] << "/data.****" << endl;

    // write the pointer to the cell variables
    for (int i = 0; i < ens_cdata_names.size(); i++)
	caseout << " scalar per element: 1  " << setw(19)
		<< setiosflags(ios::left)
		<< ens_cdata_names[i] << setw(4) << " "
		<< "./" << ens_cdata_names[i] << "/data.****" << endl;

    caseout << endl;
    // write out the time block
    caseout << "TIME" << endl;
    caseout << "time set:              " << setw(4) << "   1" << endl;
    caseout << "number of steps:       " << setw(4) << setiosflags(ios::right) 
	    << dump_times.size() << endl;
    caseout << "filename start number: " << setw(4) << "   1" << endl;
    caseout << "filename increment:    " << setw(4) << "   1" << endl;
    caseout << "time values:           " << endl;
    
    // write out times
    caseout.precision(5);
    caseout.setf(ios::scientific, ios::floatfield);
    for (int i = 0; i < dump_times.size(); i++)
	caseout << setw(12) << setiosflags(ios::right) << dump_times[i] 
		<< endl;
}

} // end of rtt_viz

//---------------------------------------------------------------------------//
//                              end of Ensight_Translator.cc
//---------------------------------------------------------------------------//
