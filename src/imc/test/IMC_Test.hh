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
#include "mc/OS_Mesh.hh"
#include "mc/OS_Builder.hh"
#include "mc/General_Topology.hh"
#include "mc/RZWedge_Mesh.hh"
#include "cdi_analytic/Analytic_Models.hh"
#include "cdi_analytic/Analytic_Gray_Opacity.hh"
#include "cdi_analytic/Analytic_Multigroup_Opacity.hh"
#include "cdi_analytic/Analytic_EoS.hh"
#include "cdi/CDI.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <iostream>
#include <vector>
#include <string>
#include <utility>

namespace rtt_mc
{
template<class PT> class Particle_Containers;
}

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

template<class PT>
class IMC_Flat_Interface :
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
    // sp to OS_Builder
    rtt_dsxx::SP<rtt_mc::OS_Builder> builder;

    // data for the Opacity and Mat_State
    SP_Data    mat_data;
    sf_double  density;
    sf_double  temperature;
    double     implicitness;
    double     delta_t;

    // data for the source builder
    double    elapsed_t;
    sf_double evol_ext;
    sf_double rad_source;
    sf_double rad_temp;
    sf_double ss_temp;
    sf_string ss_desc;

  public:
    // constructor
    IMC_Flat_Interface(rtt_dsxx::SP<rtt_mc::OS_Builder>, bool = false);

    // general interface
    double get_delta_t() const { return delta_t; }
    int    get_hybrid_diffusion_method() const { return 0; }
    
    // public interface for Opacity_Builder
    SP_Data   get_flat_data_container() const { return mat_data; }
    sf_double get_density() const {return density;}
    sf_double get_temperature() const {return temperature;}
    double    get_implicitness_factor() const { return implicitness; }

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
// IMC_FLAT_INTERFACE DEFINITIONS
//---------------------------------------------------------------------------//

// constructor
template<class PT>
IMC_Flat_Interface<PT>::IMC_Flat_Interface(
    rtt_dsxx::SP<rtt_mc::OS_Builder> osb, 
    bool                             common_mg_opacities) 
    : builder(osb),
      mat_data(new rtt_imc::Flat_Data_Container),
      density(6), 
      temperature(6),
      implicitness(1.0), 
      delta_t(.001),
      elapsed_t(.001),
      evol_ext(6),
      rad_source(6),
      rad_temp(6),
      ss_temp(2),
      ss_desc(2, "standard")
{   
    // make the Opacity and Mat_State stuff
    
    // size the data in the flat data container
    mat_data->gray_absorption_opacity.resize(6);
    mat_data->gray_scattering_opacity.resize(6);
    mat_data->mg_absorption_opacity.resize(6, sf_double(3));
    mat_data->mg_scattering_opacity.resize(6, sf_double(3));
    mat_data->specific_heat.resize(6);
    mat_data->group_boundaries.resize(4);

    // make group boundaries
    mat_data->group_boundaries[0] = 0.01;
    mat_data->group_boundaries[1] = 0.1;
    mat_data->group_boundaries[2] = 15.0;
    mat_data->group_boundaries[3] = 100.0;

    for (int i = 0; i < 3; i++)
    {
	// density
	density[i]   = 1.0;
	density[i+3] = 2.0;

	// absorption opacity in /cm
	mat_data->gray_absorption_opacity[i]    = .1  * density[i];
	mat_data->gray_absorption_opacity[i+3]  = .01 * density[i+3];

	if (common_mg_opacities)
	{
	    mat_data->mg_absorption_opacity[i][0]   = 0.1 * density[i];
	    mat_data->mg_absorption_opacity[i][1]   = 0.1 * density[i];
	    mat_data->mg_absorption_opacity[i][2]   = 0.1 * density[i];
	    
	    mat_data->mg_absorption_opacity[i+3][0] = 0.01 * density[i+3];
	    mat_data->mg_absorption_opacity[i+3][1] = 0.01 * density[i+3];
	    mat_data->mg_absorption_opacity[i+3][2] = 0.01 * density[i+3];

	    mat_data->mg_scattering_opacity[i][0]   = 0.5 * density[i];
	    mat_data->mg_scattering_opacity[i][1]   = 0.5 * density[i];
	    mat_data->mg_scattering_opacity[i][2]   = 0.5 * density[i];
	    mat_data->mg_scattering_opacity[i+3][0] = 0.0;
	    mat_data->mg_scattering_opacity[i+3][1] = 0.0;
	    mat_data->mg_scattering_opacity[i+3][2] = 0.0;
	}
	else
	{
	    mat_data->mg_absorption_opacity[i][0]   = 1.0;
	    mat_data->mg_absorption_opacity[i][1]   = 0.5;
	    mat_data->mg_absorption_opacity[i][2]   = 0.1;
	    
	    mat_data->mg_absorption_opacity[i+3][0] = 2.0;
	    mat_data->mg_absorption_opacity[i+3][1] = 1.5;
	    mat_data->mg_absorption_opacity[i+3][2] = 1.1;

	    mat_data->mg_scattering_opacity[i][0]   = 0.0;
	    mat_data->mg_scattering_opacity[i][1]   = 0.0;
	    mat_data->mg_scattering_opacity[i][2]   = 0.0;
	    mat_data->mg_scattering_opacity[i+3][0] = 0.0;
	    mat_data->mg_scattering_opacity[i+3][1] = 0.0;
	    mat_data->mg_scattering_opacity[i+3][2] = 0.0;
	}
	
	// scattering opacity in /cm
	mat_data->gray_scattering_opacity[i]    = 0.5 * density[i];
	mat_data->gray_scattering_opacity[i+3]  = 0.0 * density[i+3];

	// specific heat in jks/g/keV
	mat_data->specific_heat[i]   = .1;
	mat_data->specific_heat[i+3] = .2;

	// temperature
	temperature[i]   = 10.0;
	temperature[i+3] = 20.0;
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

//---------------------------------------------------------------------------//

template<class PT>
std::vector<std::vector<int> > IMC_Flat_Interface<PT>::get_defined_surcells()
    const
{
    return builder->get_defined_surcells();
}

//---------------------------------------------------------------------------//
//make a CDI interface for a 6 cell mesh; the data specifications for each
//cell are in Vol III pg. 1; this interface is used to test
//CDI_Mat_State_Builder only.  It contains both multigroup and gray data in
//the CDIs. The source parts of the interface have no meaning
//---------------------------------------------------------------------------//

template<class PT>
class IMC_CDI_Interface :
	public rtt_imc::Interface<PT>,
	public rtt_imc::CDI_Data_Interface
{
  public:
    typedef rtt_dsxx::SP<rtt_cdi::CDI>                       SP_CDI;
    typedef std::vector<int>                                 sf_int;
    typedef std::vector<SP_CDI>                              sf_CDI;
    typedef std::pair<int, int>                              model_pair;
    typedef std::vector<model_pair>                          sf_model_pair;
    typedef std::vector<double>                              sf_double;
    typedef std::string                                      std_string;
    typedef std::vector<std_string>                          sf_string;
    typedef std::vector<sf_int>                              vf_int;
    typedef typename rtt_mc::Particle_Containers<PT>::Census Census;
    typedef rtt_dsxx::SP<Census>                             SP_Census;


  private:
    // data for the Opacity and Mat_State
    sf_double  density;
    sf_double  temperature;
    double     implicitness;
    double     delta_t;

    sf_int        cdi_map;
    sf_CDI        cdi_list;
    sf_model_pair cdi_models;

  public:
    // constructor -> the default processor capacity is 6 cells
    IMC_CDI_Interface();

    // general interface
    double get_delta_t() const { return delta_t; }
    int    get_hybrid_diffusion_method() const { return 0; }
    
    // public interface for Opacity_Builder
    sf_double get_density() const {return density;}
    sf_double get_temperature() const {return temperature;}
    double    get_implicitness_factor() const { return implicitness; }

    // CDI specific parts of interface
    sf_CDI get_CDIs() const { return cdi_list; }
    sf_int get_CDI_map() const { return cdi_map; }
    sf_model_pair get_CDI_models() const { return cdi_models; }

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

//---------------------------------------------------------------------------//
// IMC_CDI_INTERFACE DEFINITIONS
//---------------------------------------------------------------------------//

// constructor
template<class PT>
IMC_CDI_Interface<PT>::IMC_CDI_Interface() 
    : 
      density(6),   
      temperature(6, 3.0), 
      implicitness(1.0), 
      delta_t(.001),
      cdi_map(6),
      cdi_list(3),
      cdi_models(3)
{  
    using rtt_cdi_analytic::Analytic_Gray_Opacity;
    using rtt_cdi_analytic::Analytic_Multigroup_Opacity;
    using rtt_cdi_analytic::Analytic_EoS;
    using rtt_cdi_analytic::Analytic_Opacity_Model;
    using rtt_cdi_analytic::Analytic_EoS_Model;
    using rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model;
    using rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model;
    using rtt_cdi_analytic::Constant_Analytic_Opacity_Model;
    using rtt_cdi::CDI;
    using rtt_cdi::GrayOpacity;
    using rtt_cdi::EoS;
    using rtt_cdi::MultigroupOpacity;
    using rtt_dsxx::SP;

    // make material data
    density[0] = 1.0;
    density[1] = 2.0;
    density[2] = 1.0;
    density[3] = 3.0;
    density[4] = 1.0;
    density[5] = 2.0;

    // make cdi map (see Vol III pg 1, TME)
    cdi_map[0] = 1;
    cdi_map[1] = 3;
    cdi_map[2] = 1;
    cdi_map[3] = 2;
    cdi_map[4] = 1;
    cdi_map[5] = 2;

    // make 3 cdi materials
    //  GRAY SPECIFICATION: 
    //      MAT 1: sigma = 100/T^3 cm^2/g, Cv = .1e6    kJ/g/keV
    //      MAT 2: sigma = 1.5     cm^2/g, Cv = .2e6    kJ/g/keV
    //      MAT 3: sigma = 1+.1T   cm^2/g, Cv = .1e6T^3 kJ/g/keV
    //  MULTIGROUP SPECIFICATION:
    //      Groups: = .01 .1 1 10 keV (3 groups)
    //      MAT 1: sigma = 5.0/1.5/0.5 cm^2/g
    //      MAT 2: sigma = 4.0/1.4/0.4 cm^2/g
    //      MAT 3: sigma = 3.0/1.3/0.3 cm^2/g
    // Note: 1kJ = 1e-6 Jerks

    // gray analytic opacity models
    SP<Analytic_Opacity_Model> model_1(new Polynomial_Analytic_Opacity_Model
				       (0.0, 100.0, -3.0, 0.0));
    SP<Analytic_Opacity_Model> model_2(new Constant_Analytic_Opacity_Model
				       (1.5));
    SP<Analytic_Opacity_Model> model_3(new Polynomial_Analytic_Opacity_Model
				       (1.0, 0.1, 1.0, 0.0));

    SP<Analytic_Opacity_Model> model_s(new Constant_Analytic_Opacity_Model
				       (0.0));

    // multigroup analytic opacity models
    std::vector<SP<Analytic_Opacity_Model> > mg_1(3);
    std::vector<SP<Analytic_Opacity_Model> > mg_2(3);
    std::vector<SP<Analytic_Opacity_Model> > mg_3(3);
    std::vector<SP<Analytic_Opacity_Model> > mg_s(3);

    mg_1[0] = new Constant_Analytic_Opacity_Model(5.0);
    mg_1[1] = new Constant_Analytic_Opacity_Model(1.5);
    mg_1[2] = new Constant_Analytic_Opacity_Model(0.5);

    mg_2[0] = new Constant_Analytic_Opacity_Model(4.0);
    mg_2[1] = new Constant_Analytic_Opacity_Model(1.4);
    mg_2[2] = new Constant_Analytic_Opacity_Model(0.4);

    mg_3[0] = new Constant_Analytic_Opacity_Model(3.0);
    mg_3[1] = new Constant_Analytic_Opacity_Model(1.3);
    mg_3[2] = new Constant_Analytic_Opacity_Model(0.3);

    mg_s[0] = model_s;
    mg_s[1] = model_s;
    mg_s[2] = model_s;

    std::vector<double> groups(4);
    groups[0] = 0.01;
    groups[1] = 0.1;
    groups[2] = 1.0;
    groups[3] = 10.0;

    // analytic eos models
    SP<Analytic_EoS_Model> aeos_1(
	new Polynomial_Specific_Heat_Analytic_EoS_Model(0.1e6, 0.0, 0.0, 
							0.0, 0.0, 0.0));
    SP<Analytic_EoS_Model> aeos_2(
	new Polynomial_Specific_Heat_Analytic_EoS_Model(0.2e6, 0.0, 0.0, 
							0.0, 0.0, 0.0));
    SP<Analytic_EoS_Model> aeos_3(
	new Polynomial_Specific_Heat_Analytic_EoS_Model(0.0, 0.1e6, 3.0, 
							0.0, 0.0, 0.0));

    // make gray opacities
    SP<const GrayOpacity> gop_1(new Analytic_Gray_Opacity(
				    model_1, rtt_cdi::ABSORPTION));
    SP<const GrayOpacity> gop_2(new Analytic_Gray_Opacity(
				    model_2, rtt_cdi::ABSORPTION));
    SP<const GrayOpacity> gop_3(new Analytic_Gray_Opacity(
				    model_3, rtt_cdi::ABSORPTION));
    SP<const GrayOpacity> gop_s(new Analytic_Gray_Opacity(
				    model_s, rtt_cdi::SCATTERING));

    // make multigroup opacities
    SP<const MultigroupOpacity> mgop_1(new Analytic_Multigroup_Opacity(
					   groups, mg_1, rtt_cdi::ABSORPTION));
    SP<const MultigroupOpacity> mgop_2(new Analytic_Multigroup_Opacity(
					   groups, mg_2, rtt_cdi::ABSORPTION));
    SP<const MultigroupOpacity> mgop_3(new Analytic_Multigroup_Opacity(
					   groups, mg_3, rtt_cdi::ABSORPTION));
    SP<const MultigroupOpacity> mgop_s(new Analytic_Multigroup_Opacity(
					   groups, mg_s, rtt_cdi::SCATTERING));


    // make EoS
    SP<const EoS> eos_1(new Analytic_EoS(aeos_1));
    SP<const EoS> eos_2(new Analytic_EoS(aeos_2));
    SP<const EoS> eos_3(new Analytic_EoS(aeos_3));

    // Make CDIs
    cdi_list[0] = new CDI;
    cdi_list[1] = new CDI;
    cdi_list[2] = new CDI;

    cdi_list[0]->setGrayOpacity(gop_1);
    cdi_list[0]->setGrayOpacity(gop_s);
    cdi_list[0]->setMultigroupOpacity(mgop_1);
    cdi_list[0]->setMultigroupOpacity(mgop_s);
    cdi_list[0]->setEoS(eos_1);

    cdi_list[1]->setGrayOpacity(gop_2);
    cdi_list[1]->setGrayOpacity(gop_s);
    cdi_list[1]->setMultigroupOpacity(mgop_2);
    cdi_list[1]->setMultigroupOpacity(mgop_s);
    cdi_list[1]->setEoS(eos_2);

    cdi_list[2]->setGrayOpacity(gop_3);
    cdi_list[2]->setGrayOpacity(gop_s);
    cdi_list[2]->setMultigroupOpacity(mgop_3);
    cdi_list[2]->setMultigroupOpacity(mgop_s);
    cdi_list[2]->setEoS(eos_3);

    // make model list; NOTE: here we are cheating a little because we are
    // dumping both gray and multigroup data into the same CDI's; we aren't
    // in trouble because the multigroup and gray opacities are the same
    // model; in real calculations we don't use the same CDI for both gray
    // and multigroup data; even if we do, we only check one or the other in
    // a given problem; IN SUMMARY: we only need models for either gray OR
    // multigroup; for this test we get away with both!

    // check that gray and multigroup have the same models
    if (gop_1->getModelType() != mgop_1->getModelType()) 
	throw(rtt_dsxx::assertion("Death to this test, incompatible models!"));
    else if (gop_2->getModelType() != mgop_2->getModelType()) 
	throw(rtt_dsxx::assertion("Death to this test, incompatible models!"));
    else if (gop_3->getModelType() != mgop_3->getModelType()) 
	throw(rtt_dsxx::assertion("Death to this test, incompatible models!"));
    else if (gop_s->getModelType() != mgop_s->getModelType()) 
	throw(rtt_dsxx::assertion("Death to this test, incompatible models!"));
	      
    // now make models
    cdi_models[0].first  = gop_1->getModelType();
    cdi_models[0].second = gop_s->getModelType();

    cdi_models[1].first  = gop_2->getModelType();
    cdi_models[1].second = gop_s->getModelType();

    cdi_models[2].first  = gop_3->getModelType();
    cdi_models[2].second = gop_s->getModelType();
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
