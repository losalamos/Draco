//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/tstTransporter.hh
 * \author Thomas M. Evans
 * \date   Thu May 22 16:33:34 2003
 * \brief  Flat interface for Transporter testing.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_test_tstTransporter_hh
#define rtt_imc_test_tstTransporter_hh

#include "../Interface.hh"
#include "../Flat_Data_Interface.hh"
#include "../Flat_Data_Container.hh"
#include "../Hybrid_Diffusion.hh"
#include "mc/Particle_Stack.hh"
#include "mc/OS_Builder.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_imc_test
{

//===========================================================================//
// Flat interface for Random walk testing
//===========================================================================//

template<class PT>
class RW_Interface :
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

    // hybrid diffusion scheme
    int hmodel;

  public:
    // constructor
    RW_Interface(rtt_dsxx::SP<rtt_mc::OS_Builder>);

    // general interface
    double get_delta_t() const { return delta_t; }
    int    get_hybrid_diffusion_method() const { return hmodel; }
    
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
RW_Interface<PT>::RW_Interface(rtt_dsxx::SP<rtt_mc::OS_Builder> osb) 
    : builder(osb),
      mat_data(new rtt_imc::Flat_Data_Container),
      density(osb->num_cells()), 
      temperature(osb->num_cells()),
      implicitness(1.0), 
      delta_t(.001),
      elapsed_t(.001),
      evol_ext(osb->num_cells()),
      rad_source(osb->num_cells()),
      rad_temp(osb->num_cells()),
      ss_temp(1),
      ss_desc(1, "standard"),
      hmodel(rtt_imc::Hybrid_Diffusion::RANDOM_WALK)
{   
    if (osb->num_cells() % 2 != 0)
	throw rtt_dsxx::assertion("Non-even celled test mesh.");

    // num cells (on each processor)
    int nc  = osb->num_cells();
    int mod = nc / 2;

    // make the Opacity and Mat_State stuff
    
    // size the data in the flat data container
    mat_data->gray_absorption_opacity.resize(nc);
    mat_data->gray_scattering_opacity.resize(nc);
    mat_data->rosseland_opacity.resize(nc);
    mat_data->mg_absorption_opacity.resize(nc, sf_double(3));
    mat_data->mg_scattering_opacity.resize(nc, sf_double(3));
    mat_data->specific_heat.resize(nc);
    mat_data->group_boundaries.resize(4);

    // make group boundaries
    mat_data->group_boundaries[0] = 0.01;
    mat_data->group_boundaries[1] = 0.1;
    mat_data->group_boundaries[2] = 15.0;
    mat_data->group_boundaries[3] = 100.0;

    for (int i = 0; i < mod; i++)
    {
	// density
	density[i]     = 1.0;
	density[i+mod] = 2.0;

	// absorption opacity in /cm
	mat_data->gray_absorption_opacity[i]     = 100.0  * density[i];
	mat_data->gray_absorption_opacity[i+mod] = 50.0 * density[i+mod];

	mat_data->mg_absorption_opacity[i][0]     = 0.1 * density[i];
	mat_data->mg_absorption_opacity[i][1]     = 0.1 * density[i];
	mat_data->mg_absorption_opacity[i][2]     = 0.1 * density[i];
	mat_data->mg_absorption_opacity[i+mod][0] = 0.01 * density[i+mod];
	mat_data->mg_absorption_opacity[i+mod][1] = 0.01 * density[i+mod];
	mat_data->mg_absorption_opacity[i+mod][2] = 0.01 * density[i+mod];

	mat_data->mg_scattering_opacity[i][0]     = 0.5 * density[i];
	mat_data->mg_scattering_opacity[i][1]     = 0.5 * density[i];
	mat_data->mg_scattering_opacity[i][2]     = 0.5 * density[i];
	mat_data->mg_scattering_opacity[i+mod][0] = 0.0;
	mat_data->mg_scattering_opacity[i+mod][1] = 0.0;
	mat_data->mg_scattering_opacity[i+mod][2] = 0.0;
		
	// scattering opacity in /cm
	mat_data->gray_scattering_opacity[i]      = 0.5 * density[i];
	mat_data->gray_scattering_opacity[i+mod]  = 0.0 * density[i+mod];

	// specific heat in jks/g/keV
	mat_data->specific_heat[i]     = .1;
	mat_data->specific_heat[i+mod] = .2;

	// rosseland opacities
	mat_data->rosseland_opacity[i] =
	    mat_data->gray_absorption_opacity[i] +
	    mat_data->gray_scattering_opacity[i];
	mat_data->rosseland_opacity[i+mod] = 
	    mat_data->gray_absorption_opacity[i+mod] +
	    mat_data->gray_scattering_opacity[i+mod];

	// temperature
	temperature[i]   = 1.0;
	temperature[i+mod] = 2.0;
    }

    // make the Source_Builder stuff

    if (osb->get_ss_pos().size() != 1)
	throw rtt_dsxx::assertion("Need to have 1 surface source.");

    for (int i = 0; i < nc; i++)
    {
	evol_ext[i]   = 1.0;
	rad_source[i] = 1.0;
	rad_temp[i]   = 0.0;
    }

    ss_temp[0] = 1.0;
}

//---------------------------------------------------------------------------//

template<class PT>
std::vector<std::vector<int> > RW_Interface<PT>::get_defined_surcells() const
{
    return builder->get_defined_surcells();
}

} // end namespace rtt_imc_test

#endif // rtt_imc_test_tstTransporter_hh

//---------------------------------------------------------------------------//
//              end of imc/test/tstTransporter.hh
//---------------------------------------------------------------------------//
