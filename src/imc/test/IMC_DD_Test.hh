//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/IMC_DD_Test.hh
 * \author Todd J. Urbatsch
 * \date   Tue May 30 12:49:35 2000
 * \brief  mesh, topology, interface for 9-cell DD mesh on 4 processors.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __test_IMC_DD_Test_hh__
#define __test_IMC_DD_Test_hh__

#include "../Interface.hh"
#include "../Particle.hh"
#include "mc/Coord_sys.hh"
#include "mc/XYCoord_sys.hh"
#include "mc/Layout.hh"
#include "mc/OS_Mesh.hh"
#include "mc/General_Topology.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>
#include <string>

namespace rtt_imc_dd_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

inline bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

//===========================================================================//
// BUILD DD MESHES
//===========================================================================//
// build 9-cell mesh:  2 cells on proc 0; 2 cells on proc 1; 2 cells on proc
// 2; 3 cells on proc 3.  Copied directly from rtt_mc_test.

inline rtt_dsxx::SP<rtt_mc::OS_Mesh> build_DD_Mesh() 
{
    using rtt_mc::OS_Mesh;
    using rtt_mc::XYCoord_sys;
    using rtt_mc::Layout;
    using rtt_mc::Coord_sys;
    using rtt_dsxx::SP;

    Require (C4::nodes() == 4);

    SP<OS_Mesh> mesh;

    if (C4::node() == 0)
    {
	SP<Coord_sys> cs(new XYCoord_sys());

	// build layout
	Layout lay(2);
	for (int i = 1; i <= lay.num_cells(); i++)
	    lay.set_size(i, 4);

	lay(1, 1) = 0;
	lay(1, 2) = 2;
	lay(1, 3) = 0;
	lay(1, 4) = -2;

	lay(2, 1) = 1;
	lay(2, 2) = -1;
	lay(2, 3) = 0;
	lay(2, 4) = -3;

	// build coordinate data
	OS_Mesh::vf_double vertex(2);
	vertex[0].resize(6);
	vertex[1].resize(6);
	
	vertex[0][0] = 0;
	vertex[0][1] = 1;
	vertex[0][2] = 2;
	vertex[0][3] = 0;
	vertex[0][4] = 1;
	vertex[0][5] = 2;
	vertex[1][0] = 0;
	vertex[1][1] = 0;
	vertex[1][2] = 0;
	vertex[1][3] = 1;
	vertex[1][4] = 1;
	vertex[1][5] = 1;

	OS_Mesh::vf_int cell_pair(2);
	cell_pair[0].resize(4);
	cell_pair[1].resize(4);
	
	cell_pair[0][0] = 1;
	cell_pair[0][1] = 2;
	cell_pair[0][2] = 4;
	cell_pair[0][3] = 5;
	cell_pair[1][0] = 2;
	cell_pair[1][1] = 3;
	cell_pair[1][2] = 5;
	cell_pair[1][3] = 6;
	
	mesh = new OS_Mesh(cs, lay, vertex, cell_pair, true);
    }
    else if (C4::node() == 1)
    {
	SP<Coord_sys> cs(new XYCoord_sys());

	// build layout
	Layout lay(2);
	for (int i = 1; i <= lay.num_cells(); i++)
	    lay.set_size(i, 4);

	lay(1, 1) = -2;
	lay(1, 2) = 0;
	lay(1, 3) = 0;
	lay(1, 4) = -4;

	lay(2, 1) = 0;
	lay(2, 2) = -3;
	lay(2, 3) = -1;
	lay(2, 4) = -5;

	// build coordinate data
	OS_Mesh::vf_double vertex(2);
	vertex[0].resize(8);
	vertex[1].resize(8);
	
	vertex[0][0] = 2;
	vertex[0][1] = 3;
	vertex[0][2] = 0;
	vertex[0][3] = 1;
	vertex[0][4] = 2;
	vertex[0][5] = 3;
	vertex[0][6] = 0;
	vertex[0][7] = 1;
	vertex[1][0] = 0;
	vertex[1][1] = 0;
	vertex[1][2] = 1;
	vertex[1][3] = 1;
	vertex[1][4] = 1;
	vertex[1][5] = 1;
	vertex[1][6] = 2;
	vertex[1][7] = 2;

	OS_Mesh::vf_int cell_pair(2);
	cell_pair[0].resize(4);
	cell_pair[1].resize(4);
	
	cell_pair[0][0] = 1;
	cell_pair[0][1] = 2;
	cell_pair[0][2] = 5;
	cell_pair[0][3] = 6;
	cell_pair[1][0] = 3;
	cell_pair[1][1] = 4;
	cell_pair[1][2] = 7;
	cell_pair[1][3] = 8;
	
	mesh = new OS_Mesh(cs, lay, vertex, cell_pair, true);
    }
    else if (C4::node() == 2)
    {
	SP<Coord_sys> cs(new XYCoord_sys());

	// build layout
	Layout lay(2);
	for (int i = 1; i <= lay.num_cells(); i++)
	    lay.set_size(i, 4);

	lay(1, 1) = -3;
	lay(1, 2) = 2;
	lay(1, 3) = -1;
	lay(1, 4) = -4;

	lay(2, 1) = 1;
	lay(2, 2) = 0;
	lay(2, 3) = -2;
	lay(2, 4) = -5;

	// build coordinate data
	OS_Mesh::vf_double vertex(2);
	vertex[0].resize(6);
	vertex[1].resize(6);
	
	vertex[0][0] = 1;
	vertex[0][1] = 2;
	vertex[0][2] = 3;
	vertex[0][3] = 1;
	vertex[0][4] = 2;
	vertex[0][5] = 3;
	vertex[1][0] = 1;
	vertex[1][1] = 1;
	vertex[1][2] = 1;
	vertex[1][3] = 2;
	vertex[1][4] = 2;
	vertex[1][5] = 2;

	OS_Mesh::vf_int cell_pair(2);
	cell_pair[0].resize(4);
	cell_pair[1].resize(4);
	
	cell_pair[0][0] = 1;
	cell_pair[0][1] = 2;
	cell_pair[0][2] = 4;
	cell_pair[0][3] = 5;
	cell_pair[1][0] = 2;
	cell_pair[1][1] = 3;
	cell_pair[1][2] = 5;
	cell_pair[1][3] = 6;
	
	mesh = new OS_Mesh(cs, lay, vertex, cell_pair, true);
    }
    else if (C4::node() == 3)
    {
	SP<Coord_sys> cs(new XYCoord_sys());

	// build layout
	Layout lay(3);
	for (int i = 1; i <= lay.num_cells(); i++)
	    lay.set_size(i, 4);

	lay(1, 1) = 0;
	lay(1, 2) = 2;
	lay(1, 3) = -1;
	lay(1, 4) = -0;

	lay(2, 1) = 1;
	lay(2, 2) = 3;
	lay(2, 3) = -2;
	lay(2, 4) = 0;

	lay(3, 1) = 2;
	lay(3, 2) = 0;
	lay(3, 3) = -3;
	lay(3, 4) = 0;

	// build coordinate data
	OS_Mesh::vf_double vertex(2);
	vertex[0].resize(8);
	vertex[1].resize(8);
	
	vertex[0][0] = 0;
	vertex[0][1] = 1;
	vertex[0][2] = 2;
	vertex[0][3] = 3;
	vertex[0][4] = 0;
	vertex[0][5] = 1;
	vertex[0][6] = 2;
	vertex[0][7] = 3;
	vertex[1][0] = 2;
	vertex[1][1] = 2;
	vertex[1][2] = 2;
	vertex[1][3] = 2;
	vertex[1][4] = 3;
	vertex[1][5] = 3;
	vertex[1][6] = 3;
	vertex[1][7] = 3;

	OS_Mesh::vf_int cell_pair(3);
	cell_pair[0].resize(4);
	cell_pair[1].resize(4);
	cell_pair[2].resize(4);
	
	cell_pair[0][0] = 1;
	cell_pair[0][1] = 2;
	cell_pair[0][2] = 5;
	cell_pair[0][3] = 6;
	cell_pair[1][0] = 2;
	cell_pair[1][1] = 3;
	cell_pair[1][2] = 6;
	cell_pair[1][3] = 7;
	cell_pair[2][0] = 3;
	cell_pair[2][1] = 4;
	cell_pair[2][2] = 7;
	cell_pair[2][3] = 8;
	
	mesh = new OS_Mesh(cs, lay, vertex, cell_pair, true);
    }

    return mesh;
}

//===========================================================================//
// BUILD DD TOPOLOGIES
//===========================================================================//
// build a General Topology for 9 cell mesh described above.
// copied directly from rtt_mc_test.

inline rtt_dsxx::SP<rtt_mc::Topology> build_DD_Topology()
{
    using rtt_mc::Topology;
    using rtt_mc::General_Topology;
    using rtt_dsxx::SP;
    using std::vector;

    Require (C4::nodes() == 4);

    Topology::vf_int cpp(4);
    Topology::vf_int ppc(9, vector<int>(1));
    Topology::vf_int bc(4);

    // cells-per-proc data
    cpp[0].resize(2);
    cpp[0][0] = 1;
    cpp[0][1] = 2;

    cpp[1].resize(2);
    cpp[1][0] = 3;
    cpp[1][1] = 4;

    cpp[2].resize(2);
    cpp[2][0] = 5;
    cpp[2][1] = 6;

    cpp[3].resize(3);
    cpp[3][0] = 7;
    cpp[3][1] = 8;
    cpp[3][2] = 9;

    // procs_per_cell data
    ppc[0][0] = 0;
    ppc[1][0] = 0;
    ppc[2][0] = 1;
    ppc[3][0] = 1;
    ppc[4][0] = 2;
    ppc[5][0] = 2;
    ppc[6][0] = 3;
    ppc[7][0] = 3;
    ppc[8][0] = 3;

    // boundary cell data
    bc[0].resize(3);
    bc[0][0] = 3;
    bc[0][1] = 4;
    bc[0][2] = 5;

    bc[1].resize(5);
    bc[1][0] = 1;
    bc[1][1] = 2;
    bc[1][2] = 5;
    bc[1][3] = 6;
    bc[1][4] = 7;

    bc[2].resize(5);
    bc[2][0] = 2;
    bc[2][1] = 3;
    bc[2][2] = 4;
    bc[2][3] = 8;
    bc[2][4] = 9;
    
    bc[3].resize(3);
    bc[3][0] = 4;
    bc[3][1] = 5;
    bc[3][2] = 6;

    SP<Topology> topology(new General_Topology(cpp, ppc, bc, "DD"));
    return topology;
}

//---------------------------------------------------------------------------//
// INTERFACE CLASS
//---------------------------------------------------------------------------//
// interface for the 9-cell DD mesh on 4 processors.

class IMC_DD_Interface :
    public rtt_imc::Interface<rtt_imc::Particle<rtt_mc::OS_Mesh> >
{
  private:

    // data for the Opacity and Mat_State
    sf_double  density;
    sf_double  kappa;
    sf_double  kappa_thomson;
    sf_double  temperature;
    sf_double  specific_heat;
    double     implicitness;
    double     delta_t;
    std_string analytic_opacity;
    std_string analytic_sp_heat;

    // data for topology
    int capacity;

    // data for the source builder
    double elapsed_t;
    sf_double evol_ext;
    sf_double rad_source;
    sf_double rad_temp;
    sf_double ss_temp;
    sf_string ss_desc;

  public:
    // constructor -> the default processor capacity is 2 cells
    inline IMC_DD_Interface(int = 2);
    
    // public interface for Opacity_Builder
    sf_double get_density() const {return density;}
    sf_double get_kappa() const {return kappa;}
    sf_double get_kappa_thomson() const {return kappa_thomson;}
    sf_double get_specific_heat() const {return specific_heat;}
    sf_double get_temperature() const {return temperature;}
    inline std::string get_analytic_opacity() const;
    inline std::string get_analytic_sp_heat() const;
    double get_implicit() const { return implicitness; }
    double get_delta_t() const { return delta_t; }

    // public interface for Topology
    int get_capacity() const { return capacity; }

    // public interface for Source_Builder
    double get_elapsed_t() const { return elapsed_t; }
    sf_double get_evol_ext() const { return evol_ext; }
    double get_rad_s_tend() const { return double(.1); }
    sf_double get_rad_source() const { return rad_source; }
    sf_double get_rad_temp() const { return rad_temp; }
    inline sf_string get_ss_pos() const;
    sf_double get_ss_temp() const { return ss_temp; }
    inline vf_int get_defined_surcells() const;
    int get_npnom() const { return int(1000); }
    int get_npmax() const { return int(1000); }
    double get_dnpdt() const { return double(0); }
    int get_cycle() const { return int(1); }
    std_string get_ss_dist() const { return "cosine"; }
    sf_string get_ss_desc() const { return ss_desc; }
    SP_Census get_census() const { return SP_Census(); }
    double get_ecen(int cell) const { return double(137.2); }
    // get_ecentot returns the global value
    double get_ecentot() const { return double(9*137.2); }
};

// constructor
IMC_DD_Interface::IMC_DD_Interface(int capacity_) 
    :  rtt_imc::Interface<rtt_imc::Particle<rtt_mc::OS_Mesh> >(),
       density(capacity_), 
       kappa(capacity_), 
       kappa_thomson(capacity_), 
       temperature(capacity_),
       specific_heat(capacity_), 
       implicitness(1.0), 
       delta_t(.001),
       analytic_opacity("straight"), 
       analytic_sp_heat("straight"),
       capacity(capacity_),
       elapsed_t(.001),
       evol_ext(capacity_),
       rad_source(capacity_),
       rad_temp(capacity_),
       ss_temp(1),
       ss_desc(1, "standard")
{   

    // check the hardwired numbers of cells for each processor
    Check (C4::nodes() == 4);
    Check ((C4::node() != 3) ? (capacity == 2) : true);
    Check ((C4::node() == 3) ? (capacity == 3) : true);

    // make the Opacity and Mat_State stuff
    for (int i = 0; i < capacity; i++)
    {
	// density
	density[i] = C4::node() + i + 1.0;

	// kappa (in cm^2/g)
	kappa[i] = 2 * C4::node() + i + 1.0;
	
	// kappa thomson
	kappa_thomson[i] = 2.0 * (C4::node() + i + 1.0);

	// specific heat
	specific_heat[i] = 3.0 * C4::node() + i + 1.0;

	// temperature
	temperature[i] = 3.0 * (C4::node() + i + 1.0);
    }

    // make the Source_Builder stuff

    for (int i = 0; i < capacity; i++)
    {
	evol_ext[i]   = 100;
	rad_source[i] = 200;
	rad_temp[i]   = 10.0;
    }

    ss_temp[0] = 20.0;
}

std::string IMC_DD_Interface::get_analytic_opacity() const 
{ 
    return analytic_opacity;
}

std::string IMC_DD_Interface::get_analytic_sp_heat() const 
{ 
    return analytic_sp_heat; 
}

std::vector<std::vector<int> > IMC_DD_Interface::get_defined_surcells() const
{
    // each processor has one surface source and one SS cell (cell 2)
    vf_int surcells(1);
    surcells[0].resize(1);
    surcells[0][0] = 2;

    return surcells;
}

std::vector<std::string> IMC_DD_Interface::get_ss_pos() const
{
    // each processor has one surface source
    sf_string ss_pos(1);

    if (C4::node() == 0)
	ss_pos[0] = "loy";
    else if (C4::node() == 1)
	ss_pos[0] = "lox";
    else if (C4::node() == 2)
	ss_pos[0] = "hix";
    else if (C4::node() == 3)
	ss_pos[0] = "hiy";
    else
	Insist (0, "Incorrect processor number!");

    return ss_pos;
}

} // end namespace rtt_imc_dd_test

#endif                          // __test_IMC_DD_Test_hh__

//---------------------------------------------------------------------------//
//                              end of test/IMC_DD_Test.hh
//---------------------------------------------------------------------------//
