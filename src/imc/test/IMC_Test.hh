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
#include "mc/OS_Mesh.hh"
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
// INTERFACE CLASS
//===========================================================================//
// make an interface for a 6 cell mesh

class IMC_Interface : public rtt_imc::Interface<rtt_mc::OS_Mesh>
{
  private:
    // Mesh data
    std_string coord;
    vf_double  fine_edge;
    sf_string  bnd;
    
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
    sf_string ss_pos;
    sf_double ss_temp;
    vf_int surcells;
    sf_string ss_desc;

  public:
    // constructor -> the default processor capacity is 6 cells
    inline IMC_Interface(int = 6);

    // public interface for Mesh
    std_string get_coordinates() const {return coord;}
    vf_double get_fine_edge() const {return fine_edge;} 
    sf_string get_boundaries() const {return bnd;}
    
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

    // public interface for Source_Init
    void set_defined_surcells(int i, const sf_int &sc) {surcells[i-1] = sc;}

    // public interface for Source_Builder
    double get_elapsed_t() const { return elapsed_t; }
    sf_double get_evol_ext() const { return evol_ext; }
    double get_rad_s_tend() const { return double(.1); }
    sf_double get_rad_source() const { return rad_source; }
    sf_double get_rad_temp() const { return rad_temp; }
    sf_string get_ss_pos() const { return ss_pos; }
    sf_double get_ss_temp() const { return ss_temp; }
    vf_int get_defined_surcells() const { return surcells; }
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
IMC_Interface::IMC_Interface(int capacity_)
    :  rtt_imc::Interface<rtt_mc::OS_Mesh>(), 
       coord("xy"), 
       fine_edge(2),
       bnd(4),
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
       ss_pos(2),
       ss_temp(2),
       surcells(2),
       ss_desc(2, "standard")
{
    // make the Mesh stuff

    // calculate the fine edges
    fine_edge[0].resize(4);
    fine_edge[1].resize(3);

    fine_edge[0][0] = -1;
    fine_edge[1][0] = -1;

    for (int i = 1; i < fine_edge[0].size(); i++)
	fine_edge[0][i] = fine_edge[0][i-1] + 1;
    for (int i = 1; i < fine_edge[1].size(); i++)
	fine_edge[1][i] = fine_edge[1][i-1] + 2;

    // calculate the boundaries
    for (int i = 1; i < 4; i++)
	bnd[i] = "vacuum";
    bnd[0] = "reflect";
    
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

    ss_pos[0]  = "loy";
    ss_pos[1]  = "hix";
    ss_temp[0] = 20.0;
    ss_temp[1] = 0.0;
    
    surcells[0].resize(2);
    surcells[0][0] = 1;
    surcells[0][1] = 2;
}


std::string IMC_Interface::get_analytic_opacity() const 
{ 
    return analytic_opacity;
}

std::string IMC_Interface::get_analytic_sp_heat() const 
{ 
    return analytic_sp_heat; 
}

} // end namespace rtt_imc_test

#endif                          // __imc_test_IMC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of imc/test/IMC_Test.hh
//---------------------------------------------------------------------------//
