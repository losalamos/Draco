//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/CDI_Mat_State_Builder.hh
 * \author Thomas M. Evans
 * \date   Fri Nov 16 11:23:17 2001
 * \brief  CDI_Mat_State_Builder class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_CDI_Mat_State_Builder_hh__
#define __imc_CDI_Mat_State_Builder_hh__

#include "Mat_State_Builder.hh"
#include "cdi/CDI.hh"
#include "ds++/Assert.hh"
#include <vector>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class CDI_Mat_State_Builder
 *
 * This class builds instances of rtt_imc::Opacity and rtt_imc::Mat_State.
 * It receives data necessary to build CDI objects--from the rtt_cdi
 * package--that read data from files or build analytic data.  It receives an
 * Interface Type (IT) in its constructor.  The IT must be a derived class of
 * rtt_imc::Interface and rtt_imc::CDI_Data_Interface.
 *
 * This builder uses the rtt_cdi package to build cell-centered material and
 * opacity data.  As such, it requires material IDs, model descriptions, and
 * the like to build the necessary rtt_cdi::CDI objects.  The member
 * functions that provide data to this builder are defined in the
 * rtt_imc::CDI_Data_Interface class.
 *
 * This class, along with rtt_imc::Flat_Mat_State_Builder, is a derived class
 * of Mat_State_Builder.  It should be used when a client is using CDI
 * (rtt_cdi) to define material data
 *
 * \sa rtt_imc::Interface, rtt_imc::CDI_Data_Interface.  See
 * tstMat_Data_Builder for examples of usage.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class CDI_Mat_State_Builder : public Mat_State_Builder<MT>
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>             SP_Mesh;
    typedef rtt_dsxx::SP<Mat_State<MT> > SP_Mat_State;
    typedef rtt_dsxx::SP<Opacity<MT> >   SP_Opacity;
    typedef rtt_dsxx::SP<rtt_cdi::CDI>   SP_CDI;
    typedef std::vector<SP_CDI>          sf_CDI;
    typedef std::vector<double>          sf_double;
    typedef std::vector<int>             sf_int;

  private:
    // Material CDI objects.
    sf_CDI    material_cdi;
    
    // Map of material_cdi to cells
    sf_int    cdi_cell_map;

    // Cell-centered densities in g/cc.
    sf_double density;
    
    // Cell-centered temperatures in keV.
    sf_double temperature;
    
    // Fleck and Cummings implicitness factor.
    double    implicitness;

    // Timestep in shakes.
    double    delta_t;

  public:
    // Constructor.
    template<class IT>
    explicit CDI_Mat_State_Builder(rtt_dsxx::SP<IT>);

    // >>> PUBLIC INTERFACE

    // Build the Mat_State.
    SP_Mat_State build_Mat_State(SP_Mesh) const;

    // Build the Opacity.
    SP_Opacity build_Opacity(SP_Mesh m, SP_Mat_State s) const;
};

//---------------------------------------------------------------------------//
// MEMBER TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
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
template<class MT>
template<class IT>
CDI_Mat_State_Builder<MT>::CDI_Mat_State_Builder(rtt_dsxx::SP<IT> interface)
    : Mat_State_Builder<MT>()
{
    Require (interface);

    // assign data members from the interface parser
    material_cdi = interface->get_CDIs();
    cdi_cell_map = interface->get_CDI_map();
    density      = interface->get_density();
    temperature  = interface->get_temperature();
    implicitness = interface->get_implicitness_factor();
    delta_t      = interface->get_delta_t();

    Ensure (delta_t > 0.0);
    Ensure (implicitness >= 0.0 && implicitness <= 1.0);
    Ensure (density.size() == temperature.size());
    Ensure (density.size() == cdi_cell_map.size());
}

} // end namespace rtt_imc

#endif                          // __imc_CDI_Mat_State_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/CDI_Mat_State_Builder.hh
//---------------------------------------------------------------------------//
