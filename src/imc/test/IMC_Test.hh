//----------------------------------*-C++-*----------------------------------//
// IMC_Test.hh
// Thomas M. Evans
// Tue Apr 27 10:53:05 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Some services we will need to test the imc packages
//---------------------------------------------------------------------------//

#ifndef __imc_test_IMC_Test_hh__
#define __imc_test_IMC_Test_hh__

#include "../Interface.hh"
#include "../Particle.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/General_Topology.hh"
#include "mc/RZWedge_Mesh.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>
#include <string>

namespace rtt_imc_test
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
// PARSER CLASS FOR MESH BUILDERS
//===========================================================================//
// make a simple parser class that tells mesh builders the name of the mesh
// input file

class Parser
{
  public:
    std::string file_name;
    explicit Parser(std::string fn) : file_name(fn) {/*...*/}
    std::string get_mesh_file() const { return file_name; }
};

//===========================================================================//
// INTERFACE CLASS
//===========================================================================//
// make an interface for a 6 cell mesh

class IMC_Interface : 
    public rtt_imc::Interface<rtt_imc::Particle<rtt_mc::OS_Mesh> >
{
  private:
    // sp to OS_Builder
    rtt_dsxx::SP<rtt_mc::OS_Builder> builder;

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
    // constructor -> the default processor capacity is 6 cells
    inline IMC_Interface(rtt_dsxx::SP<rtt_mc::OS_Builder>, int = 6);
    
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
    sf_string get_ss_pos() const { return builder->get_ss_pos(); }
    sf_double get_ss_temp() const { return ss_temp; }
    inline vf_int get_defined_surcells() const;
    int get_npnom() const { return int(1000); }
    int get_npmax() const { return int(1000); }
    double get_dnpdt() const { return double(0); }
    int get_cycle() const { return int(1); }
    std_string get_ss_dist() const { return "cosine"; }
    sf_string get_ss_desc() const { return ss_desc; }
    SP_Census get_census() const { return SP_Census(); }
    double get_ecen(int cell) const { return double(0); }
    double get_ecentot() const { return double(0); }
};

// constructor
IMC_Interface::IMC_Interface(rtt_dsxx::SP<rtt_mc::OS_Builder> osb, 
			     int capacity_) 
    :  rtt_imc::Interface<rtt_imc::Particle<rtt_mc::OS_Mesh> >(),
       builder(osb),
       density(6), 
       kappa(6), 
       kappa_thomson(6), 
       temperature(6),
       specific_heat(6), 
       implicitness(1.0), 
       delta_t(.001),
       analytic_opacity("straight"), 
       analytic_sp_heat("straight"),
       capacity(capacity_),
       elapsed_t(.001),
       evol_ext(6),
       rad_source(6),
       rad_temp(6),
       ss_temp(2),
       ss_desc(2, "standard")
{   
    // make the Opacity and Mat_State stuff

    for (int i = 0; i < 3; i++)
    {
	// density
	density[i]   = 1.0;
	density[i+3] = 2.0;

	// kappa (in cm^2/g)
	kappa[i]     = .1;
	kappa[i+3]   = .01;
	
	// kappa thomson
	kappa_thomson[i]   = .5;
	kappa_thomson[i+3] = 0.0;

	// specific heat
	specific_heat[i]   = .1;
	specific_heat[i+3] = .2;

	// temperature
	temperature[i]   = 10;
	temperature[i+3] = 20;
    }

    // make the Source_Builder stuff

    for (int i = 0; i < 6; i++)
    {
	evol_ext[i]   = 100;
	rad_source[i] = 200;
	rad_temp[i]   = 10.0;
    }

    ss_temp[0] = 20.0;
    ss_temp[1] = 0.0;
}


std::string IMC_Interface::get_analytic_opacity() const 
{ 
    return analytic_opacity;
}

std::string IMC_Interface::get_analytic_sp_heat() const 
{ 
    return analytic_sp_heat; 
}

std::vector<std::vector<int> > IMC_Interface::get_defined_surcells() const
{
    return builder->get_defined_surcells();
}

//===========================================================================//
// MAKE AN AMR RZWEDGE_MESH
//===========================================================================//

rtt_dsxx::SP<rtt_mc::RZWedge_Mesh> make_RZWedge_Mesh_AMR(double);

//===========================================================================//
// SEND/RECEIVE GENERAL DD TOPOLOGY
//===========================================================================//
// send out a DD Topology

inline void send_TOP(const rtt_mc::General_Topology &topology)
{
    Require (C4::node() == 0);

    // get General Topology data
    rtt_mc::Topology::vf_int cpp = topology.get_cells_per_proc();
    rtt_mc::Topology::vf_int ppc = topology.get_procs_per_cell();
    rtt_mc::Topology::vf_int bc  = topology.get_bound_cells();

    int size = 0;
    for (int i = 0; i < cpp.size(); i++)
	size += cpp[i].size();
    for (int i = 0; i < ppc.size(); i++)
	size += ppc[i].size();
    for (int i = 0; i < bc.size(); i++)
	size += bc[i].size();
	    
    int *idata = new int[size];
    int *jdata = new int[cpp.size() + ppc.size() + bc.size()];
	
    int ic = 0;
    int jc = 0;
    for (int i = 0; i < cpp.size(); i++)
    {
	jdata[jc++] = cpp[i].size();
	for (int j = 0; j < cpp[i].size(); j++)
	    idata[ic++] = cpp[i][j];
    }

    for (int i = 0; i < ppc.size(); i++)
    {
	jdata[jc++] = ppc[i].size();
	for (int j = 0; j < ppc[i].size(); j++)
	    idata[ic++] = ppc[i][j];
    }

    for (int i = 0; i < bc.size(); i++)
    {
	jdata[jc++] = bc[i].size();
	for (int j = 0; j < bc[i].size(); j++)
	    idata[ic++] = bc[i][j];
    }

    Check (ic == size);
    Check (jc == cpp.size() + ppc.size() + bc.size());

    C4::Send<int>(cpp.size(), 1, 100);
    C4::Send<int>(ppc.size(), 1, 101);
    C4::Send<int>(bc.size(),  1, 102);
    C4::Send<int>(size,       1, 103);
    C4::Send<int>(jdata, jc, 1, 104);
    C4::Send<int>(idata, ic, 1, 105);

    delete [] idata;
    delete [] jdata;   
}

//---------------------------------------------------------------------------//
// receive a DD_Topology

inline rtt_dsxx::SP<rtt_mc::Topology> recv_TOP() 
{
    Require (C4::node());

    int cpp_size;
    int ppc_size;
    int bc_size;
    int size;

    C4::Recv(cpp_size, 0, 100);
    C4::Recv(ppc_size, 0, 101);
    C4::Recv(bc_size,  0, 102);
    C4::Recv(size,     0, 103);

    int *idata = new int[size];
    int *jdata = new int[cpp_size + ppc_size + bc_size];

    C4::Recv(jdata, cpp_size + ppc_size + bc_size, 0, 104);
    C4::Recv(idata, size, 0, 105);

    rtt_mc::Topology::vf_int cpp(cpp_size);
    rtt_mc::Topology::vf_int ppc(ppc_size);
    rtt_mc::Topology::vf_int bc(bc_size);

    int jc = 0; 
    int ic = 0;
    for (int i = 0; i < cpp.size(); i++)
    {
	cpp[i].resize(jdata[jc++]);
	for (int j = 0; j < cpp[i].size(); j++)
	    cpp[i][j] = idata[ic++];
    }

    for (int i = 0; i < ppc.size(); i++)
    {
	ppc[i].resize(jdata[jc++]);
	for (int j = 0; j < ppc[i].size(); j++)
	    ppc[i][j] = idata[ic++];
    }

    for (int i = 0; i < bc.size(); i++)
    {
	bc[i].resize(jdata[jc++]);
	for (int j = 0; j < bc[i].size(); j++)
	    bc[i][j] = idata[ic++];
    }

    Check (ic == size);
    Check (jc == cpp_size + ppc_size + bc_size);

    rtt_dsxx::SP<rtt_mc::Topology> return_top;
    return_top = new rtt_mc::General_Topology(cpp,ppc,bc,"DD");

    delete [] idata;
    delete [] jdata;

    return return_top;
}

} // end namespace rtt_imc_test

#endif                          // __imc_test_IMC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of imc/test/IMC_Test.hh
//---------------------------------------------------------------------------//
