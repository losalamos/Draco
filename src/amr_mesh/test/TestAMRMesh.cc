//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   amr_mesh/test/TestAMRMesh.cc
 * \author B.T. Adams
 * \date   Tue Mar 14 09:48:00 2000
 * \brief  Header file for the RTT_Format class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestAMRMesh.hh"
#include "../Release.hh"
#include "UnitTestFrame/PassFailStream.hh"
#include <sstream>
#include <iostream>
#include <vector>

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_amr_test::TestAMRMesh;
    
    return SP<TestApp>(new TestAMRMesh(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_amr_test
{

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::multimap;
using std::make_pair;
using std::set;
using std::pair;
using std::ostream;
using std::fill;
using rtt_dsxx::SP;
using rtt_amr::CAR_CU_Interface;
using rtt_meshReaders::RTT_Format;

TestAMRMesh::TestAMRMesh(int argc, char * argv[], ostream & os_in)
    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestAMRMesh" << endl;
}

string TestAMRMesh::version() const
{
    return rtt_amr::release();
}

/*!
 * \brief Tests the RTT_format mesh reader.
 *
 */
string TestAMRMesh::runTest()
{
    // New meshes added to this test will have to be added to the enumeration
    // Meshes in the header file.
    const int MAX_MESHES = 1;
    string filename[MAX_MESHES] = {"str1"};
    Meshes mesh_type;
    bool verbose = true;

    for (int mesh_number = 0; mesh_number < MAX_MESHES; mesh_number++)
    {
	// Construct a new CAR_CU_Interface class object.
	SP<CAR_CU_Interface> interface(new 
			   CAR_CU_Interface(filename[mesh_number], verbose));
	pass(" Construct ") << "Read " << filename[mesh_number]
			    << " and created Interface class without coreing"
			    << " in or firing an assertion." 
			    << endl;

	// Construct a new RTT_Format class object and parse the input files.
	SP<RTT_Format> rttMesh = interface->parser();
	pass(" Construct ") << "Read " << filename[mesh_number] << ".mesh"
			    << " and created RTT_Format class without coreing"
			    << " in or firing an assertion." 
			    << endl;

	bool all_passed = true;
	// The following switch allows addition of other meshes for testing,
	// with the "DEFINED" mesh providing an example. Only the check_dims
	// tests is required and it will be automatically called by the other
	// tests (with the exception of check header) if not invoked herein.
	// The comparison data must also be provided for additional meshes
	// within the switch structure residing in the test functions. 
        switch (mesh_number)
	{
	// Test all nested class accessor functions for a very simplistic 
	// mesh file (enum STR1).
	case (0):
	    mesh_type = STR1;
	    all_passed = all_passed && check_interface(* interface, * rttMesh, 
						       mesh_type);
	    break;

	default:
	    fail("runTest") << "Invalid mesh type encountered." << endl;
	    all_passed = false;
	    break;
	}
	if (!all_passed)
	    fail(filename[mesh_number]) << "Errors occured testing mesh " 
					<< "number " << mesh_type << endl;
    }

    // Report results of test.
    if (passed())
    {
	return "All tests passed.";
    }
    return "Some tests failed.";
}

bool TestAMRMesh::check_interface(const rtt_amr::CAR_CU_Interface & interface,
				  const rtt_meshReaders::RTT_Format & mesh,
				  const Meshes & meshtype)
{
    // Exercise the interface accessor functions for this mesh.
    bool all_passed = true;
    string coord_system;
    vector<double> density;
    vector<double> kappa;
    vector<double> kappa_thomson;
    vector<double> specific_heat;
    vector<double> temperature;
    string analytic_opacity;
    string analytic_sp_heat;
    double implicit;
    vector<double> rad_temp;
    vector<double> evol_ext;
    vector<double> rad_source;
    double rad_s_tend;
    vector<string> ss_pos;
    vector<double> ss_temp;
    vector<vector<int> > ss_cells;
    double delta_t;
    int npmax, npnom, dnpdt, capacity;
    string ss_dist;
    int max_cycle, printf, buffer, seed;
    string surface_file, mesh_file;

    switch (meshtype)
    {
    case STR1:
        coord_system = "XYZ";
	density.resize(mesh.get_dims_ncells(), 1.0);
	kappa.resize(mesh.get_dims_ncells(), 1.0);
	kappa_thomson.resize(mesh.get_dims_ncells(), 0.5);
	specific_heat.resize(mesh.get_dims_ncells(), 0.1);
	temperature.resize(mesh.get_dims_ncells(), 1.0);
	analytic_opacity = "straight";
	analytic_sp_heat = "straight";
	implicit = 1.0;
	rad_temp.resize(mesh.get_dims_ncells(), 1.5);
	evol_ext.resize(0);
	rad_source.resize(0);
	rad_s_tend = 0;
	ss_pos.resize(1,"lox");
	ss_temp.resize(1,1.0);
	ss_cells.resize(1);
	ss_cells[0].push_back(1);  ss_cells[0].push_back(5);
	ss_cells[0].push_back(9);  ss_cells[0].push_back(13);
	ss_cells[0].push_back(17); ss_cells[0].push_back(21);
	ss_cells[0].push_back(24); ss_cells[0].push_back(26);
	ss_cells[0].push_back(30); ss_cells[0].push_back(34);
	ss_cells[0].push_back(36); ss_cells[0].push_back(38);
	ss_cells[0].push_back(42); ss_cells[0].push_back(46);
	ss_cells[0].push_back(50); ss_cells[0].push_back(54);
	delta_t = 0.01;
	npmax = 100;
	npnom = 100;
	dnpdt = 2;
	capacity = 64;
	ss_dist = "normal";
	max_cycle = 5;
	printf = 1;
	buffer = 1;
	seed = 9347593;
	surface_file = "undefined";
	mesh_file = "str1.mesh";
	break;

   default:
        fail("check_interface") << "Invalid mesh type encountered." << endl;
	all_passed = false;
	return all_passed;
    }
    // Check coordinate system.
    if (interface.get_coordinates() != coord_system)
    {
        fail(" Interface Coordinates ") << 
	     "Interface coord_system not obtained." << endl;
 	all_passed = false;
    }
    // Check density.
    if (interface.get_density() != density)
    {
        fail(" Interface Density ") << "Interface density not obtained." 
				    << endl;
 	all_passed = false;
    }
    // Check kappa.
    if (interface.get_kappa() != kappa)
    {
        fail(" Interface Kappa ") << "Interface kappa not obtained." 
				  << endl;
 	all_passed = false;
    }
    // Check kappa thomson.
    if (interface.get_kappa_thomson() != kappa_thomson)
    {
        fail(" Interface Kappa Thomson ") << 
	     "Interface kappa_thomson not obtained." << endl;
 	all_passed = false;
    }
    // Check specific heat.
    if (interface.get_specific_heat() != specific_heat)
    {
        fail(" Interface Specific Heat ") << 
	     "Interface specific_heat not obtained." << endl;
 	all_passed = false;
    }
    // Check temperature.
    if (interface.get_temperature() != temperature)
    {
        fail(" Interface Temperature ") 
	     << "Interface temperature not obtained." << endl;
 	all_passed = false;
    }
    // Check analytic opacity.
    if (interface.get_analytic_opacity() != analytic_opacity)
    {
        fail(" Interface Analytic Opacity ") << 
	     "Interface analytic opacity not obtained." << endl;
 	all_passed = false;
    }
    // Check analytic specific heat.
    if (interface.get_analytic_sp_heat() != analytic_sp_heat)
    {
        fail(" Interface Analytic Specific Heat ") << 
	     "Interface analytic specific heat not obtained." << endl;
 	all_passed = false;
    }
    // Check implicitness factor.
    if (interface.get_implicit() != implicit)
    {
        fail(" Interface Implicitness ") << 
	     "Interface implicitness not obtained." << endl;
 	all_passed = false;
    }
    // Check zone size.
    if (interface.get_zone_size() != mesh.get_dims_ncells())
    {
        fail(" Interface Zone Size ") << 
	     "Interface zone_size not obtained." << endl;
 	all_passed = false;
    }
    // Check radiation temperature.
    if (interface.get_rad_temp() != rad_temp)
    {
        fail(" Interface Radiation Temperature ") <<
	     "Interface rad_temp not obtained." << endl;
 	all_passed = false;
    }
    // Check external volumetric source.
    if (interface.get_evol_ext() != evol_ext)
    {
        fail(" Interface Volumetric Source ") << 
	     "Interface evol_ext not obtained." << endl;
 	all_passed = false;
    }
    // Check cell external volumetric source.
    bool got_evol_ext = true;
    for (int i = 0; i < evol_ext.size(); i++)
        if (evol_ext[i] != interface.get_evol_ext(i + 1))
	    got_evol_ext = false;
    if (!got_evol_ext)
    {
        fail(" Interface Volumetric Source ") << 
	     "Interface cell evol_ext not obtained." << endl;
 	all_passed = false;
    }
    // Check external radiation source.
    if (interface.get_rad_source() != rad_source)
    {
        fail(" Interface Radiation Source ") 
	     << "Interface rad_source not obtained." << endl;
 	all_passed = false;
    }
    // Check cell external radiation source.
    bool got_rad_source = true;
    for (int i = 0; i < rad_source.size(); i++)
        if (rad_source[i] != interface.get_rad_source(i + 1))
	    got_rad_source = false;
    if (!got_rad_source)
    {
        fail(" Interface Radiation Source ") << 
	     "Interface cell rad_source not obtained." << endl;
 	all_passed = false;
    }
    // Check surface source positions size.
    if (interface.get_ss_pos_size() != ss_pos.size())
    {
        fail(" Interface Surface Source Positions ") << 
	     "Interface ss_pos_size not obtained." << endl;
 	all_passed = false;
    }
    // Check surface source positions.
    if (interface.get_ss_pos() != ss_pos)
    {
        fail(" Interface Surface Source Positions ") << 
	     "Interface ss_pos not obtained." << endl;
 	all_passed = false;
    }
    // Check surface source position.
    bool got_ss_pos = true;
    for (int i = 0; i < ss_pos.size(); i++)
        if (ss_pos[i] != interface.get_ss_pos(i + 1))
	    got_ss_pos = false;
    if (!got_ss_pos)
    {
        fail(" Interface Surface Source Positions ") << 
	     "Interface surface ss_pos not obtained." << endl;
 	all_passed = false;
    }    
     // Check surface source temperatures.
    if (interface.get_ss_temp() != ss_temp)
    {
        fail(" Interface Surface Source Temperatures ") << 
	     "Interface ss_temp not obtained." << endl;
 	all_passed = false;
    }
    // Check surface source temperature.
    bool got_ss_temp = true;
    for (int i = 0; i < ss_temp.size(); i++)
        if (ss_temp[i] != interface.get_ss_temp(i + 1))
	    got_ss_temp = false;
    if (!got_ss_temp)
    {
        fail(" Interface Surface Source Temperatures ") << 
	     "Interface surface ss_temp not obtained." << endl;
 	all_passed = false;
    }    
    // Check surface source cells size.
    bool got_ss_cells_size = true;
    for (int i = 0; i < ss_cells.size(); i++)
        if (ss_cells[i].size() != interface.get_ss_cells_size(i + 1))
	    got_ss_cells_size = false;
    if (!got_ss_cells_size)
    {
        fail(" Interface Surface Source Cells ") << 
	     "Interface ss_cells_size not obtained." << endl;
 	all_passed = false;
    }
    // Check surface source cells.
    if (ss_cells != interface.get_defined_surcells())
    {
        fail(" Interface Surface Source Cells ") << 
	     "Interface defined_surcells not obtained." << endl;
 	all_passed = false;
    }
    // Check surface source cells.
    bool got_ss_cells = true;
    for (int i = 0; i < ss_cells.size(); i++)
        if (ss_cells[i] != interface.get_defined_surcells(i + 1))
	    got_ss_cells = false;
    if (!got_ss_cells)
    {
        fail(" Interface Surface Source Cells ") << 
	     "Interface surface defined_surcells not obtained." << endl;
 	all_passed = false;
    }
    // Check external radiation source cutoff time.
    if (interface.get_rad_s_tend() != rad_s_tend)
    {
        fail(" Interface Radiation Source Cutoff Time ") << 
	     "Interface rad_s_tend not obtained." << endl;
 	all_passed = false;
    }
   // Check time step.
    if (interface.get_delta_t() != delta_t)
    {
        fail(" Interface Time Step ") << "Interface delta_t not obtained." 
				      << endl;
 	all_passed = false;
    }
     // Check maximum number of particles.
    if (interface.get_npmax() != npmax)
    {
        fail(" Interface Maximum Number of Particles ") << 
	     "Interface npmax not obtained." << endl;
 	all_passed = false;
    }
   // Check nominal number of particles.
    if (interface.get_npnom() != npnom)
    {
        fail(" Interface Nominal Number of Particles ") << 
	     "Interface npnom not obtained." << endl;
 	all_passed = false;
    }
    // Check particle time rate of change.
    if (interface.get_dnpdt() != dnpdt)
    {
        fail(" Interface Particle Time Rate of Change ") << 
	     "Interface dnpdt not obtained." << endl;
 	all_passed = false;
    }
    // Check capacity.
    if (interface.get_capacity() != capacity)
    {
        fail(" Interface Capacity ") << "Interface capacity not obtained." 
				     << endl;
 	all_passed = false;
    }
     // Check surface source distribution.
    if (interface.get_ss_dist() != ss_dist)
    {
        fail(" Interface SS Distribution ") << 
	     "Interface ss_dist not obtained." << endl;
 	all_passed = false;
    }
    // Check maximum cycles.
    if (interface.get_max_cycle() != max_cycle)
    {
        fail(" Interface Maximum Cycles ") << 
	     "Interface max_cycle not obtained." << endl;
 	all_passed = false;
    }
    // Check print frequency.
    if (interface.get_printf() != printf)
    {
        fail(" Interface Print Frequency ") << 
	     "Interface printf not obtained." << endl;
 	all_passed = false;
    }
    // Check buffer.
    if (interface.get_buffer() != buffer)
    {
        fail(" Interface Buffer ") << "Interface buffer not obtained." << endl;
 	all_passed = false;
    }
    // Check seed.
    if (interface.get_seed() != seed)
    {
        fail(" Interface Seed ") << "Interface seed not obtained." << endl;
 	all_passed = false;
    }
   // Check surface file.
    if (interface.get_surface_file() != surface_file)
    {
        fail(" Interface Surface File ") << 
	     "Interface surface file not obtained." << endl;
 	all_passed = false;
    }
    // Check mesh file.
    if (interface.get_mesh_file() != mesh_file)
    {
        fail(" Interface Mesh File ") << "Interface mesh file not obtained." 
				      << endl;
 	all_passed = false;
    }

    // Check that all Interface class accessors passed their tests.
    if (all_passed)
        pass(" Interface Accessors " ) << "Got all Interface accessors." 
				       << endl;
    else
	fail(" Interface Accessors ") << "Errors in some Interface accessors." 
				      << endl;

    return all_passed;
}

} // end namespace rtt_amr_test


//---------------------------------------------------------------------------//
//                              end of amr_mesh/test/TestAMRMesh.cc
//---------------------------------------------------------------------------//
