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
#include "mc/Topology.hh"
#include "mc/OS_Mesh.hh"
#include "ds++/SP.hh"
#include <iostream>

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

rtt_dsxx::SP<rtt_mc::OS_Mesh> build_DD_Mesh();

//===========================================================================//
// BUILD DD TOPOLOGIES
//===========================================================================//
// build a General Topology for 9 cell mesh described above.
// copied directly from rtt_mc_test.

rtt_dsxx::SP<rtt_mc::Topology> build_DD_Topology();

//===========================================================================//
// DD INTERFACE CLASS 
//===========================================================================//
// interface for the 9-cell DD mesh on 4 processors.

class IMC_DD_Interface :
    public rtt_imc::Interface<rtt_imc::Particle<rtt_mc::OS_Mesh> >
{
  private:
    // data for the Opacity and Mat_State
    sf_double  density;
    sf_double  kappa;
    sf_double  kappa_offset;
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
    IMC_DD_Interface(int = 2);
    
    // public interface for Opacity_Builder
    sf_double get_density() const {return density;}
    sf_double get_kappa() const {return kappa;}
    sf_double get_kappa_offset() const {return kappa_offset;}
    sf_double get_kappa_thomson() const {return kappa_thomson;}
    sf_double get_specific_heat() const {return specific_heat;}
    sf_double get_temperature() const {return temperature;}
    std_string get_analytic_opacity() const;
    std_string get_analytic_sp_heat() const;
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
    sf_string get_ss_pos() const;
    sf_double get_ss_temp() const { return ss_temp; }
    vf_int get_defined_surcells() const;
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

} // end namespace rtt_imc_dd_test

#endif                          // __test_IMC_DD_Test_hh__

//---------------------------------------------------------------------------//
//                              end of test/IMC_DD_Test.hh
//---------------------------------------------------------------------------//
