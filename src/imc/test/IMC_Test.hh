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
#include "../Flat_Data_Interface.hh"
#include "../CDI_Data_Interface.hh"
#include "../Particle.hh"
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/General_Topology.hh"
#include "mc/RZWedge_Mesh.hh"
#include "cdi/CDI.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>
#include <string>

namespace rtt_imc_test
{

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
// INTERFACE CLASSES
//===========================================================================//
// make a flat interface for a 6 cell mesh

class IMC_Flat_Interface :
	public rtt_imc::Interface<rtt_imc::Particle<rtt_mc::OS_Mesh> >,
	public rtt_imc::Flat_Data_Interface
{
  private:
    // sp to OS_Builder
    rtt_dsxx::SP<rtt_mc::OS_Builder> builder;

    // data for the Opacity and Mat_State
    sf_double  density;
    sf_double  absorption;
    sf_double  scattering;
    sf_double  temperature;
    sf_double  specific_heat;
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
    // constructor -> the default processor capacity is 6 cells
    IMC_Flat_Interface(rtt_dsxx::SP<rtt_mc::OS_Builder>, int = 6);
    
    // public interface for Opacity_Builder
    sf_double get_density() const {return density;}
    sf_double get_absorption_opacity() const { return absorption; }
    sf_double get_scattering_opacity() const { return scattering; }
    sf_double get_specific_heat() const {return specific_heat;}
    sf_double get_temperature() const {return temperature;}
    double    get_implicitness_factor() const { return implicitness; }
    double    get_delta_t() const { return delta_t; }

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
    vf_int get_defined_surcells() const;
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

//---------------------------------------------------------------------------//
// make a CDI interface for a 6 cell mesh; the data specifications for each
// cell are in Vol III pg. 1; this interface is used to test
// CDI_Mat_State_Builder only, the source parts of the interface have no
// meaning 

class IMC_CDI_Interface :
	public rtt_imc::Interface<rtt_imc::Particle<rtt_mc::OS_Mesh> >,
	public rtt_imc::CDI_Data_Interface
{
  public:
    typedef rtt_dsxx::SP<rtt_cdi::CDI> SP_CDI;
    typedef std::vector<int>           sf_int;
    typedef std::vector<SP_CDI>        sf_CDI;

  private:
    // data for the Opacity and Mat_State
    sf_double  density;
    sf_double  temperature;
    double     implicitness;
    double     delta_t;

    sf_int     cdi_map;
    sf_CDI     cdi_list;

  public:
    // constructor -> the default processor capacity is 6 cells
    IMC_CDI_Interface();
    
    // public interface for Opacity_Builder
    sf_double get_density() const {return density;}
    sf_double get_temperature() const {return temperature;}
    double    get_implicitness_factor() const { return implicitness; }
    double    get_delta_t() const { return delta_t; }

    // CDI specific parts of interface
    sf_CDI get_CDIs() const { return cdi_list; }
    sf_int get_CDI_map() const { return cdi_map; }

    // public interface for Topology
    int get_capacity() const { return int(); }

    // public interface for Source_Builder
    double get_elapsed_t() const { return double(); }
    sf_double get_evol_ext() const { return sf_double(); }
    double get_rad_s_tend() const { return double(); }
    sf_double get_rad_source() const { return sf_double(); }
    sf_double get_rad_temp() const { return sf_double(); }
    sf_string get_ss_pos() const { return sf_string(); }
    sf_double get_ss_temp() const { return sf_double(); }
    vf_int get_defined_surcells() const { return vf_int(); }
    int get_npnom() const { return int(); }
    int get_npmax() const { return int(); }
    double get_dnpdt() const { return double(0); }
    int get_cycle() const { return int(); }
    std_string get_ss_dist() const { return std_string(); }
    sf_string get_ss_desc() const { return sf_string(); }
    SP_Census get_census() const { return SP_Census(); }
    double get_ecen(int cell) const { return double(); }
    double get_ecentot() const { return double(); }
};

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
    using rtt_mc::Topology;

    Require (C4::node() == 0);

    // Pack up the General Topology
    rtt_dsxx::SP<Topology::Pack> pack = topology.pack();

    C4::Send<int>(pack->get_parallel_scheme_indicator(), 1, 100);
    C4::Send<int>(pack->get_size(), 1, 101);
    C4::Send<int>(pack->begin(), pack->get_size(), 1, 102);
}

//---------------------------------------------------------------------------//
// receive a DD_Topology

inline rtt_dsxx::SP<rtt_mc::Topology> recv_TOP() 
{
    using rtt_mc::General_Topology;

    Require (C4::node());

    int indicator;
    int size;

    C4::Recv(indicator, 0, 100);
    C4::Recv(size, 0, 101);

    int *data = new int[size];
    C4::Recv<int>(data, size, 0, 102);

    General_Topology::Pack pack(indicator, size, data);

    rtt_dsxx::SP<rtt_mc::Topology> return_top = pack.unpack();

    return return_top;
}

} // end namespace rtt_imc_test

#endif                          // __imc_test_IMC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of imc/test/IMC_Test.hh
//---------------------------------------------------------------------------//
