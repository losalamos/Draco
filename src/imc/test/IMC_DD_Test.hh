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
#include "../Flat_Data_Interface.hh"
#include "../Particle.hh"
#include "mc/Topology.hh"
#include "mc/OS_Mesh.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>
#include <string>

namespace rtt_mc
{
template<class PT> class Particle_Containers;
}

namespace rtt_imc_dd_test
{

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

template<class PT>
class IMC_DD_Interface :
	public rtt_imc::Interface<PT>,
	public rtt_imc::Flat_Data_Interface
{
  public:
    // Useful typedefs
    typedef std::vector<double>                              sf_double;
    typedef std::vector<int>                                 sf_int;
    typedef std::vector<sf_int>                              vf_int;
    typedef std::vector<std::string>                         sf_string;
    typedef std::string                                      std_string;
    typedef rtt_dsxx::SP<rtt_imc::Flat_Data_Container>       SP_Data;
    typedef typename rtt_mc::Particle_Containers<PT>::Census Census;
    typedef rtt_dsxx::SP<Census>                             SP_Census;

  private:
    // data for the Opacity and Mat_State
    SP_Data    mat_data;
    sf_double  density;
    sf_double  temperature;
    double     implicitness;
    double     delta_t;

    // data for topology
    int capacity;

    // data for the source builder
    double    elapsed_t;
    sf_double evol_ext;
    sf_double rad_source;
    sf_double rad_temp;
    sf_double ss_temp;
    sf_string ss_desc;

  public:
    // constructor -> the default processor capacity is 2 cells
    IMC_DD_Interface(int = 2);
    
    // public interface for Opacity_Builder
    SP_Data   get_flat_data_container() const {return mat_data;}
    sf_double get_density() const {return density;}
    sf_double get_temperature() const {return temperature;}
    double    get_implicitness_factor() const { return implicitness; }
    double    get_delta_t() const { return delta_t; }

    // processor capacity
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

//---------------------------------------------------------------------------//
// DD INTERFACE CLASS MEMBER DEFINITIONS
//---------------------------------------------------------------------------//
// constructor

template<class PT>
IMC_DD_Interface<PT>::IMC_DD_Interface(int capacity_) 
    : mat_data(new rtt_imc::Flat_Data_Container), 
      density(capacity_),
      temperature(capacity_),
      implicitness(1.0), 
      delta_t(.001),
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
    
    // size the data in the flat data container
    mat_data->gray_absorption_opacity.resize(capacity);
    mat_data->gray_scattering_opacity.resize(capacity);
    mat_data->specific_heat.resize(capacity);

    for (int i = 0; i < capacity; i++)
    {
	// density
	density[i] = C4::node() + i + 1.0;

	// absorption in /cm
	mat_data->gray_absorption_opacity[i] = (2 * C4::node() + i + 1.0) *
	    density[i];

	// scattering in /cm
	mat_data->gray_scattering_opacity[i] = (2.0 * (C4::node() + i + 1.0)) * 
	    density[i];

	// specific heat
	mat_data->specific_heat[i] = 3.0 * C4::node() + i + 1.0;

	// temperature
	mat_data->temperature[i]   = 3.0 * (C4::node() + i + 1.0);
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

//---------------------------------------------------------------------------//

template<class PT>
std::vector<std::vector<int> > IMC_DD_Interface<PT>::get_defined_surcells()
    const
{
    // each processor has one surface source and one SS cell (cell 2)
    vf_int surcells(1);
    surcells[0].resize(1);
    surcells[0][0] = 2;

    return surcells;
}

//---------------------------------------------------------------------------//

template<class PT>
std::vector<std::string> IMC_DD_Interface<PT>::get_ss_pos() const
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
