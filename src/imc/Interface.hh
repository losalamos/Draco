//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   imc/Interface.hh
 * \author Thomas M. Evans
 * \date   Wed Jul  7 15:42:05 1999 
 * \brief  Interface abstract base class definition for the imc package.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Interface_hh__
#define __imc_Interface_hh__

#include "Particle_Buffer.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Interface
 *
 * \brief Interface base class for multi-cycle imc components.
 *
 * The Interface base class defines interface functionality required by the
 * Topology_Builder, Source_Builder, Flat_Mat_State_Builder,
 * CDI_Mat_State_Builder classes.  The primary use of this class is to build
 * interfaces that can provide these factories with the proper data to build
 * Opacity, Mat_State, and Source.  These are all components in the rtt_imc
 * package. Additionally, these are all components that change from timestep
 * to timestep.  For example, a MT may or may not change between cycles.
 * However, the Mat_State is always changed.  Thus. various classes may take
 * on the responsibility for becoming interfaces.  To accomplish this, simply
 * inherit the Interface base class and the necessary public interface for
 * Opacity_Builder, Topology_Builder, and Source_Builder will be enforced at
 * compile time.
 *
 * The Mat_State_Builder derived classes: Flat_Mat_State_Builder and
 * CDI_Mat_State_Builder, require different interfaces.  The functions that
 * are common to both reside within rtt_imc::Interface.  The functions that
 * are unique to Flat_Mat_State_Builder are defined in the
 * Flat_Data_Interface virtual class.  The functions that are unique to
 * CDI_Mat_State_Builder are defined in the CDI_Data_Interface virtual class.
 * To make an interface class that will support the imc package using cdi,
 * the interface must inherit (or contain the functions declarations) from
 * rtt_imc::Interface and rtt_imc::CDI_Data_Interface.  There are no name
 * collisions between rtt_imc::Interface and rtt_imc::Flat_Data_Interface and
 * rtt_cdi::CDI_Data_Interface.
 *
 * \sa rtt_imc::Flat_Data_Interface and rtt_imc::CDI_Data_Interface.
 *
 * Peruse the tests tstMat_State_Builder.cc, tstSource_Builder.cc, and
 * tstTopology_Builder.cc for examples.
 */
// revision history:
// -----------------
// 0) original
// 1) 31 Jul 2001 : Added kappa_offset (tju)
// 2) 13 Nov 2001 : modified interface class to conform to new
//                  Mat_State_Builder constructions.
// 
//===========================================================================//

template<class PT>
class Interface 
{
  public:
    //! Useful typedefs for derived interfaces.
    typedef std::string                                        std_string;
    typedef std::vector<double>                                sf_double;
    typedef std::vector<std::vector<double> >                  vf_double;
    typedef std::vector<int>                                   sf_int;
    typedef std::vector<std::vector<int> >                     vf_int;
    typedef std::vector<std::string>                           sf_string;
    typedef rtt_dsxx::SP<typename Particle_Buffer<PT>::Census> SP_Census; 
    
  public:
    // Constructor.
    Interface() { /* no data to construct */ }

    // Virtual constructor to make life happy down the inhertiance chain.
    virtual ~Interface() { /* need a destructor for inheritance chain */ }
   
    // >>> FUNCTIONS REQUIRED BY MULTIPLE BUILDER CLASSES
    
    //! Get cell-centered densities in g/cc.
    virtual sf_double get_density() const = 0;

    //! Get cell-centered temperatures in keV.
    virtual sf_double get_temperature() const = 0;

    //! Get Fleck and Cummings implicitness factor.
    virtual double    get_implicitness_factor() const = 0;

    //! Get the timestep in shakes.
    virtual double    get_delta_t() const = 0;

    // >>> FUNCTIONS REQUIRED BY TOPOLOGY BUILDER
    
    // Get the cells per processor capacity.
    virtual int       get_capacity() const = 0;
    
    // >>> FUNCTIONS REQUIRED BY SOURCE_BUILDERs

    // Interface to Source_Builder base class.

    //! Get elapsed time to beginning of timestep of imc run in shakes.
    virtual double     get_elapsed_t()        const = 0;

    //! Get cell-centered external material source in Jerks/cm^3/shake.
    virtual sf_double  get_evol_ext()         const = 0;

    //! Get turn-off time of radiation source in shakes.
    virtual double     get_rad_s_tend()       const = 0;
    
    //! Get cell-centered radiation source in Jerks/cm^3/shake.
    virtual sf_double  get_rad_source()       const = 0;

    //! Get cell-centered radiation temperature in keV.
    virtual sf_double  get_rad_temp()         const = 0;

    //! Get vector of surface source positions.
    virtual sf_string  get_ss_pos()           const = 0;

    //! Get vector of surface source temperatures in keV.
    virtual sf_double  get_ss_temp()          const = 0;

    //! Get vector field of surface source cell indices.
    virtual vf_int     get_defined_surcells() const = 0;

    //! Get nominal number of particles for this timestep.
    virtual int        get_npnom()            const = 0;

    //! Get maximum number of particles for this timestep.
    virtual int        get_npmax()            const = 0;

    //! Get differential number of particle/(shake of elapsed time).
    virtual double     get_dnpdt()            const = 0;

    //! Get problem cycle number.
    virtual int        get_cycle()            const = 0;

    //! Get the surface source angular distribution descriptor.
    virtual std_string get_ss_dist()          const = 0;

    //! Get the surface source descriptor.
    virtual sf_string  get_ss_desc()          const = 0;

    // persistent source data, this data is calculated and is persistent
    // between time steps; at the beginning of a non-initial timestep this
    // data must be retrieved from storage

    // Get the particle census.
    virtual SP_Census  get_census()           const = 0;

    // Get the census energy in a cell in keV.
    virtual double     get_ecen(int cell)     const = 0;

    // Get the total, global census energy in keV.
    virtual double     get_ecentot()          const = 0;
};

} // end namespace rtt_imc

#endif                          // __imc_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Interface.hh
//---------------------------------------------------------------------------//
