//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/IMC_Test.cc
 * \author Thomas M. Evans
 * \date   Thu Aug  3 14:20:15 2000
 * \brief  Definitions from IMC_Test.hh.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "../Global.hh"
#include "mc/AMR_Layout.hh"
#include "mc/XYZCoord_sys.hh"
#include "cdi_analytic/Analytic_Models.hh"
#include "cdi_analytic/Analytic_Gray_Opacity.hh"
#include "cdi_analytic/Analytic_EoS.hh"
#include <cmath>

namespace rtt_imc_test
{

using rtt_mc::global::pi;
using rtt_mc::RZWedge_Mesh;
using rtt_mc::AMR_Layout;
using rtt_mc::Coord_sys;
using rtt_mc::XYZCoord_sys;
using rtt_cdi_analytic::Analytic_Gray_Opacity;
using rtt_cdi_analytic::Analytic_EoS;
using rtt_cdi_analytic::Analytic_Opacity_Model;
using rtt_cdi_analytic::Analytic_EoS_Model;
using rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model;
using rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model;
using rtt_cdi_analytic::Constant_Analytic_Opacity_Model;
using rtt_cdi::CDI;
using rtt_cdi::GrayOpacity;
using rtt_cdi::EoS;
using rtt_dsxx::SP;

using std::sqrt;
using std::cos;
using std::sin;

typedef std::vector<double>    sf_double;
typedef std::vector<sf_double> vf_double;

//===========================================================================//
// INTERFACE CLASSES
//===========================================================================//

//---------------------------------------------------------------------------//
// IMC_FLAT_INTERFACE DEFINITIONS
//---------------------------------------------------------------------------//

// constructor
IMC_Flat_Interface::IMC_Flat_Interface(rtt_dsxx::SP<rtt_mc::OS_Builder> osb, 
				       int capacity_) 
    : builder(osb),
      density(6), 
      absorption(6), 
      scattering(6),  
      temperature(6),
      specific_heat(6), 
      implicitness(1.0), 
      delta_t(.001),
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

	// absorption opacity in /cm
	absorption[i]     = .1  * density[i];
	absorption[i+3]   = .01 * density[i+3];
	
	// scattering opacity in /cm
	scattering[i]   = .5  * density[i];
	scattering[i+3] = 0.0 * density[i+3];

	// specific heat in jks/g/keV
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

//---------------------------------------------------------------------------//

std::vector<std::vector<int> > IMC_Flat_Interface::get_defined_surcells()
    const
{
    return builder->get_defined_surcells();
}

//---------------------------------------------------------------------------//
// IMC_CDI_INTERFACE DEFINITIONS
//---------------------------------------------------------------------------//

// constructor
IMC_CDI_Interface::IMC_CDI_Interface() 
    : 
      density(6),   
      temperature(6, 3.0), 
      implicitness(1.0), 
      delta_t(.001),
      cdi_map(6),
      cdi_list(3)
{  
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

    // make 3 cdi materials: 
    //      MAT 1: sigma = 100/T^3 cm^2/g, Cv = .1e6    kJ/g/keV
    //      MAT 2: sigma = 1.5     cm^2/g, Cv = .2e6    kJ/g/keV
    //      MAT 3: sigma = 1+.1T   cm^2/g, Cv = .1e6T^3 kJ/g/keV
    // Note: 1kJ = 1e-6 Jerks

    // analytic opacity models
    SP<Analytic_Opacity_Model> model_1(new Polynomial_Analytic_Opacity_Model
				       (0.0, 100.0, -3.0, 1.0, 0.0));
    SP<Analytic_Opacity_Model> model_2(new Constant_Analytic_Opacity_Model
				       (1.5));
    SP<Analytic_Opacity_Model> model_3(new Polynomial_Analytic_Opacity_Model
				       (1.0, 0.1, 1.0, 1.0, 0.0));

    SP<Analytic_Opacity_Model> model_s(new Constant_Analytic_Opacity_Model
				       (0.0));

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
    cdi_list[0]->setEoS(eos_1);

    cdi_list[1]->setGrayOpacity(gop_2);
    cdi_list[1]->setGrayOpacity(gop_s);
    cdi_list[1]->setEoS(eos_2);

    cdi_list[2]->setGrayOpacity(gop_3);
    cdi_list[2]->setGrayOpacity(gop_s);
    cdi_list[2]->setEoS(eos_3);
}

//===========================================================================//
// BUILD A RZWEDGE MESH
//===========================================================================//

SP<RZWedge_Mesh> make_RZWedge_Mesh_AMR(double phi)
{
    double rphi = phi * pi / 180.0;
    
    // make coord
    SP<Coord_sys> coord(new XYZCoord_sys());

    // make layout (single level mesh)
    AMR_Layout layout(9, 6);
    {
	layout(1, 1, 1) = 1;
	layout(1, 2, 1) = 4;
	layout(1, 3, 1) = 1;
	layout(1, 4, 1) = 1;
	layout(1, 5, 1) = 1;
	layout(1, 6, 1) = 2;
    
	layout(2, 1, 1) = 2;
	layout(2, 2, 1) = 5;
	layout(2, 3, 1) = 2;
	layout(2, 4, 1) = 2;
	layout(2, 5, 1) = 1;
	layout(2, 6, 1) = 3;
    
	layout.set_size(3, 2, 2);
	layout(3, 1, 1) = 3;
	layout(3, 2, 1) = 6;
	layout(3, 2, 2) = 7;
	layout(3, 3, 1) = 3;
	layout(3, 4, 1) = 3;
	layout(3, 5, 1) = 2;
	layout(3, 6, 1) = 0;
    
	layout(4, 1, 1) = 1;
	layout(4, 2, 1) = 0;
	layout(4, 3, 1) = 4;
	layout(4, 4, 1) = 4;
	layout(4, 5, 1) = 4;
	layout(4, 6, 1) = 5;
    
	layout.set_size(5, 6, 2);
	layout(5, 1, 1) = 2;
	layout(5, 2, 1) = 0;
	layout(5, 3, 1) = 5;
	layout(5, 4, 1) = 5;
	layout(5, 5, 1) = 4;
	layout(5, 6, 1) = 6;
	layout(5, 6, 2) = 8;
    
	layout(6, 1, 1) = 3;
	layout(6, 2, 1) = 8;
	layout(6, 3, 1) = 6;
	layout(6, 4, 1) = 6;
	layout(6, 5, 1) = 5;
	layout(6, 6, 1) = 7;
    
	layout(7, 1, 1) = 3;
	layout(7, 2, 1) = 9;
	layout(7, 3, 1) = 7;
	layout(7, 4, 1) = 7;
	layout(7, 5, 1) = 6;
	layout(7, 6, 1) = 0;
    
	layout(8, 1, 1) = 6;
	layout(8, 2, 1) = 0;
	layout(8, 3, 1) = 8;
	layout(8, 4, 1) = 8;
	layout(8, 5, 1) = 5;
	layout(8, 6, 1) = 9;
    
	layout(9, 1, 1) = 7;
	layout(9, 2, 1) = 0;
	layout(9, 3, 1) = 9;
	layout(9, 4, 1) = 9;
	layout(9, 5, 1) = 8;
	layout(9, 6, 1) = 0;
    }

    // make x-z extents for each cell
    vf_double xz(9, sf_double(4));
    vf_double rz(9, sf_double(4));
    {
	double rconv = sqrt(rphi/sin(rphi))*cos(rphi/2);

	rz[0][0] = 0/rconv;
	rz[0][1] = 1/rconv;
	rz[0][2] = 0;
	rz[0][3] = 1;
	
	rz[1][0] = 0/rconv;
	rz[1][1] = 1/rconv;
	rz[1][2] = 1;
	rz[1][3] = 2;
	
	rz[2][0] = 0/rconv;
	rz[2][1] = 1/rconv;
	rz[2][2] = 2;
	rz[2][3] = 3;
	
	rz[3][0] = 1/rconv;
	rz[3][1] = 2/rconv;
	rz[3][2] = 0;
	rz[3][3] = 1;
	
	rz[4][0] = 1/rconv;
	rz[4][1] = 2/rconv;
	rz[4][2] = 1;
	rz[4][3] = 2;
	
	rz[5][0] = 1/rconv;
	rz[5][1] = 1.5/rconv;
	rz[5][2] = 2;
	rz[5][3] = 2.5;
	
	rz[6][0] = 1/rconv;
	rz[6][1] = 1.5/rconv;
	rz[6][2] = 2.5;
	rz[6][3] = 3;
	
	rz[7][0] = 1.5/rconv;
	rz[7][1] = 2/rconv;
	rz[7][2] = 2;
	rz[7][3] = 2.5;
	
	rz[8][0] = 1.5/rconv;
	rz[8][1] = 2/rconv;
	rz[8][2] = 2.5;
	rz[8][3] = 3;

	for (int i = 0; i < layout.num_cells(); i++)
	{
	    // convert rz faces
	    for (int j = 0; j < 2; j++)
	    {
		xz[i][j]   = sqrt(rphi/sin(rphi)) * rz[i][j] * cos(rphi/2.0); 
		xz[i][j+2] = rz[i][j+2];
	    }
	}
    }

    // make the mesh
    SP<RZWedge_Mesh> mesh(new RZWedge_Mesh(coord, layout, xz, phi));
    Check (mesh->full_Mesh());

    return mesh;
}

} // end namespace rtt_imc_test


//---------------------------------------------------------------------------//
//                              end of IMC_Test.cc
//---------------------------------------------------------------------------//
