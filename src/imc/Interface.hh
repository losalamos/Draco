//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   imc/Interface.hh
 * \author Thomas M. Evans
 * \date   Wed Jul  7 15:42:05 1999 
 * \brief  Interface abstract base class definition
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Interface_hh__
#define __imc_Interface_hh__

#include "Particle.hh"
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
 * The interface base class defines the interface functionality required by
 * the Opacity_Builder, Topology_Builder, and Source_Builder classes.  The
 * primary use of this class is to build interfaces that can provide these
 * factories with the proper data to build Opacity, Mat_State, and Source.
 * These are all components in the rtt_imc package.  Additionally, these are
 * all components that change from timestep to timestep.  For example, a MT
 * may or may not change between cycles.  However, the Mat_State is always
 * changed.  Thus. various classes may take on the responsibility for
 * becoming interfaces.  To accomplish this, simply inherit the Interface
 * base class and the necessary public interface for Opacity_Builder,
 * Topology_Builder, and Source_Builder will be enforced at compile time.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, class PT = Particle<MT> >
class Interface 
{
  public:
    // useful typedefs
    typedef std::string                                        std_string;
    typedef std::vector<double>                                sf_double;
    typedef std::vector<std::vector<double> >                  vf_double;
    typedef std::vector<int>                                   sf_int;
    typedef std::vector<std::vector<int> >                     vf_int;
    typedef std::vector<std::string>                           sf_string;
    typedef rtt_dsxx::SP<typename Particle_Buffer<PT>::Census> SP_Census; 
    
  public:
    // constructor
    Interface() { /* no data to construct */ }

    // virtual constructor to make life happy down the inhertiance chain
    virtual ~Interface() { /* need a destructor for inheritance chain */ }

    // >>> FUNCTIONS REQUIRED BY OPACITY_BUILDER

    // Material data.
    virtual sf_double get_density() const = 0;
    virtual sf_double get_temperature() const = 0;

    // Strings that determine how to calculate the opacity and specific
    // heats.
    virtual std_string get_analytic_opacity() const = 0;
    virtual std_string get_analytic_sp_heat() const = 0;

    // Opacity data.
    virtual sf_double get_kappa() const = 0;
    virtual sf_double get_kappa_thomson() const = 0;
    virtual sf_double get_specific_heat() const = 0;

    // Fleck implicitness factor.
    virtual double get_implicit() const = 0;
   
    // Timestep.
    virtual double get_delta_t() const = 0;

    // >>> FUNCTIONS REQUIRED BY TOPOLOGY BUILDER
    
    // Get the cells per processor capacity.
    virtual int get_capacity() const = 0;
    
    // >>> FUNCTIONS REQUIRED BY SOURCE_BUILDERs

    // Interface to Source_Builder base class.
    virtual double get_elapsed_t() const = 0;
    virtual sf_double get_evol_ext() const = 0;
    virtual double get_rad_s_tend() const = 0;
    virtual sf_double get_rad_source() const = 0;
    virtual sf_double get_rad_temp() const = 0;
    virtual sf_string get_ss_pos() const = 0;
    virtual sf_double get_ss_temp() const = 0;
    virtual vf_int get_defined_surcells() const = 0;
    virtual int get_npnom() const = 0;
    virtual int get_npmax() const = 0;
    virtual double get_dnpdt() const = 0;
    virtual int get_cycle() const = 0;
    virtual std_string get_ss_dist() const = 0;
    virtual sf_string get_ss_desc() const = 0;

    // persistent source data, this data is calculated and is persistent
    // between time steps; at the beginning of a non-initial timestep this
    // data must be retrieved from storage
    virtual SP_Census get_census() const = 0;
    virtual double get_ecen(int cell) const = 0;
    virtual double get_ecentot() const = 0;
};

} // end namespace rtt_imc

#endif                          // __imc_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Interface.hh
//---------------------------------------------------------------------------//
