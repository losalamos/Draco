//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/CDI_Mat_State_Builder.hh
 * \author Thomas M. Evans
 * \date   Fri Nov 16 11:23:17 2001
 * \brief  CDI_Mat_State_Builder class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_CDI_Mat_State_Builder_hh
#define rtt_imc_CDI_Mat_State_Builder_hh

#include "Mat_State_Builder.hh"
#include "Frequency.hh"
#include "Opacity.hh"
#include "Mat_State.hh"
#include "Diffusion_Opacity.hh"
#include "Hybrid_Diffusion.hh"
#include "Global.hh"
#include "cdi/CDI.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <vector>

namespace rtt_imc
{

// Forward declarations.
template<class MT> class Fleck_Factors;

//===========================================================================//
/*!
 * \class CDI_Mat_State_Builder
 *
 * This class builds instances of rtt_imc::Gray_Frequency,
 * rtt_imc::Multigroup_Frequency, rtt_imc::Opacity,
 * rtt_imc::Diffusion_Opacity, and rtt_imc::Mat_State.  It receives data
 * necessary to build CDI objects--from the rtt_cdi package--that read data
 * from files or build analytic data.  It receives an Interface Type (IT) in
 * its constructor.  The IT must be a derived class of rtt_imc::Interface and
 * rtt_imc::CDI_Data_Interface.
 *
 * This builder uses the rtt_cdi package to build cell-centered material and
 * opacity data.  As such, it requires material IDs, model descriptions, and
 * the like to build the necessary rtt_cdi::CDI objects.  The member
 * functions that provide data to this builder are defined in the
 * rtt_imc::CDI_Data_Interface class.
 *
 * The CDIs \b must contain data for absorption and scattering. The
 * rtt_cdi::Model types for each are specified by the interface. If
 * Diffusion_Opacity objects need to be built then the data required depends
 * upon the frequency treatment.  For gray, if the CDIs have the Rosseland
 * total opacity set (rtt_cdi::Model == rtt_cdi::ROSSELAND and
 * rtt_cdi::Reaction == rtt_cdi::TOTAL) then that is used.  Otherwise, the
 * absorption and scattering opacities are added together to estimate the
 * Rosseland total opacity.  However, if the absorption opacity model is
 * rtt_cdi::PLANCK then an assertion is thrown.
 *
 * For multigroup, the Rosseland opacities are calculated by integrating over
 * each group as follows:
 * \f[
 * \frac{1}{\sigma_{R}} = \frac{\int(\sigma_{a}+\sigma_{s})
 * \frac{\partial B}{\partial T}}{\int\frac{\partial B}{\partial T}}
 * \f]
 *
 * This class, along with rtt_imc::Flat_Mat_State_Builder, is a derived class
 * of Mat_State_Builder.  It should be used when a client is using CDI
 * (rtt_cdi) to define material data
 *
 * \sa rtt_imc::Interface, rtt_imc::CDI_Data_Interface,
 * rtt_imc::Mat_State_Builder.  See imc/test/tstCDIMat_Data_Builder.cc for
 * examples of usage.
 */
/*!
 * \example imc/test/tstCDI_Mat_State_Builder.cc
 *
 * Test of CDI_Mat_State_Builder derived type.
 */
// revision history:
// -----------------
// 0) original
// 1) 20-FEB-2003 : updated to match new Mat_State_Builder interface
// 2) 06-MAR-2003 : updated to build diffusion opacities for gray problems
// 3) 08-AUG-2003 : updated to build diffusion opacities for mg problems
// 4) 11-AUG-2003 : updated to use partial template specialization (like
//                  Flat_Mat_State_Builder) instead of class specialization.
// 
//===========================================================================//

template<class MT, class FT>
class CDI_Mat_State_Builder : public Mat_State_Builder<MT,FT>
{
  public:
    // Useful typedefs.

    // Frequency-type independent objects.
    typedef rtt_dsxx::SP<MT>                      SP_Mesh;
    typedef rtt_dsxx::SP<Mat_State<MT> >          SP_Mat_State;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >         SP_Opacity;
    typedef rtt_dsxx::SP<FT>                      SP_Frequency;
    typedef rtt_dsxx::SP<Diffusion_Opacity<MT> >  SP_Diff_Opacity;
    typedef rtt_dsxx::SP<Fleck_Factors<MT> >      SP_Fleck_Factors;
    typedef rtt_dsxx::SP<rtt_cdi::CDI>            SP_CDI;
    typedef std::vector<SP_CDI>                   sf_CDI;
    typedef std::pair<int, int>                   model_pair;
    typedef std::vector<model_pair>               sf_model_pair;
    typedef std::vector<double>                   sf_double;
    typedef std::vector<sf_double>                vf_double;
    typedef std::vector<int>                      sf_int;

    // Frequency dependent switches for partial specialization.
    typedef rtt_imc::global::Type_Switch<Gray_Frequency>       Switch_Gray;
    typedef rtt_imc::global::Type_Switch<Multigroup_Frequency> Switch_MG;
    typedef rtt_dsxx::SP<Opacity<MT,Gray_Frequency> >          SP_Gray_Opacity;
    typedef rtt_dsxx::SP<Opacity<MT,Multigroup_Frequency> >    SP_MG_Opacity;
    typedef rtt_dsxx::SP<Gray_Frequency>                       SP_Gray;
    typedef rtt_dsxx::SP<Multigroup_Frequency>                 SP_MG;
    typedef typename rtt_imc::global::Type_Switch<FT>::Type    Dummy_Type;
    typedef typename MT::template CCSF<sf_double>              ccvf;

  private:
    // Material CDI objects.
    sf_CDI        material_cdi;
    
    // Map of material_cdi to cells
    sf_int        cdi_cell_map;

    // List of (absorption,scattering) models for each CDI.
    sf_model_pair cdi_models;

    // Cell-centered densities in g/cc.
    sf_double     density;
    
    // Cell-centered temperatures in keV.
    sf_double     temperature;
    
    // Fleck and Cummings implicitness factor.
    double        implicitness;

    // Timestep in shakes.
    double        delta_t;

    // Switch for building diffusion opacities.
    bool          build_diffusion_opacity;

  private:
    // >>> BUILT OBJECTS
    
    // Built frequency.
    SP_Frequency frequency;

    // Built Mat_State.
    SP_Mat_State mat_state;

    // Built Opacity.
    SP_Opacity opacity;
    
    // Built Diffusion_Opacity.
    SP_Diff_Opacity diff_opacity;

  private:
    // >>> IMPLEMENTATION FUNCTIONS

    // Build the Mat_State.
    void build_mat_state(SP_Mesh);

  private:
    // >>> PARTIAL SPECIALIZATIONS ON FREQUENCY TYPE

    // Build a Gray_Frequency.
    template<class Stop_Explicit_Instantiation>
    void build_frequency(Switch_Gray);

    // Build a Multigroup_Frequency.
    template<class Stop_Explicit_Instantiation>
    void build_frequency(Switch_MG);

    // Build an Opacity<MT,Gray_Frequency>.
    template<class Stop_Explicit_Instantiation>
    void build_opacity(Switch_Gray, SP_Mesh);

    // Build a Diffusion_Opacity when FT=Gray_Frequency.
    template<class Stop_Explicit_Instantiation>
    void build_diff_opacity_gray(Switch_Gray, SP_Mesh, SP_Fleck_Factors);

    // Build an Opacity<MT,Multigroup_Frequency>.
    template<class Stop_Explicit_Instantiation>
    void build_opacity(Switch_MG, SP_Mesh);

    // Build MG absorption and scattering opacities.
    template<class Stop_Explicit_Instantiation>
    void build_mg_opacities(Switch_MG, ccvf &, ccvf &);

  public:
    // Constructor.
    template<class IT>
    explicit CDI_Mat_State_Builder(rtt_dsxx::SP<IT>);

    // >>> PUBLIC INTERFACE

    // Build frequency, material state, and opacity.
    void build_mat_classes(SP_Mesh);

    //! Get the frequency.
    SP_Frequency get_Frequency() const { return frequency; }
    
    //! Get the mat state.
    SP_Mat_State get_Mat_State() const { return mat_state; }

    //! Get the opacity
    SP_Opacity get_Opacity() const { return opacity; }

    //! Get the diffusion opacity.
    SP_Diff_Opacity get_Diffusion_Opacity() const { return diff_opacity; }
};

//---------------------------------------------------------------------------//
// MEMBER TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for CDI_Mat_State_Builder.
 *
 * The constructor gets data from the interface type.  The interface type
 * must provide member functions defined by rtt_imc::Interface and
 * rtt_imc::CDI_Data_Interface.  However, because this is used as a template
 * argument, it is not required that the interface type inherit from these
 * classes.  It is only required that the interface type have these functions
 * defined.
 *
 * \param interface rtt_dsxx::SP to an interface type that contains the
 * interface specification defined by rtt_imc::Interface and
 * rtt_imc::CDI_Data_Interface. 
 */
template<class MT, class FT>
template<class IT>
CDI_Mat_State_Builder<MT,FT>::CDI_Mat_State_Builder(rtt_dsxx::SP<IT> interface)
    : Mat_State_Builder<MT,FT>(),
      material_cdi(interface->get_CDIs()),
      cdi_cell_map(interface->get_CDI_map()),
      cdi_models(interface->get_CDI_models()),
      density(interface->get_density()),
      temperature(interface->get_temperature()),
      implicitness(interface->get_implicitness_factor()),
      delta_t(interface->get_delta_t()),
      build_diffusion_opacity(false)
{
    Require (interface);
    
    // set switch
    int hybrid = interface->get_hybrid_diffusion_method();
    switch (hybrid)
    {
    case Hybrid_Diffusion::TRANSPORT:
	build_diffusion_opacity = false;
	break;

    case Hybrid_Diffusion::RANDOM_WALK:
    case Hybrid_Diffusion::DDIMC:
	build_diffusion_opacity = true;
	break;

    default:
	throw rtt_dsxx::assertion("Invalid hybrid diffusion scheme.");
	break;
    }

    Ensure (delta_t > 0.0);
    Ensure (implicitness >= 0.0 && implicitness <= 1.0);
    Ensure (density.size() == temperature.size());
    Ensure (density.size() == cdi_cell_map.size());
    Ensure (material_cdi.size() > 0);
    Ensure (cdi_models.size() == material_cdi.size());
}

} // end namespace rtt_imc

#endif                          // rtt_imc_CDI_Mat_State_Builder_hh

//---------------------------------------------------------------------------//
//                              end of imc/CDI_Mat_State_Builder.hh
//---------------------------------------------------------------------------//
