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

#include "ds++/SP.hh"
#include "matprops/InterpedMaterialProps.hh"
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
 * the Opacity_Builder and Source_Init (Parallel_Source_Init) classes.  The
 * primary use of this class is to build interfaces that can provide these
 * factories with the proper data to build Opacity, Mat_State, and Source.
 * These are all components in the rtt_imc package.  Additionally, these are
 * all components that change from timestep to timestep.  For example, a MT
 * may or may not change between cycles.  However, the Mat_State is always
 * changed.  Thus. various classes may take on the responsibility for
 * becoming interfaces.  To accomplish this, simply inherit the Interface
 * base class and the necessary public interface for Opacity_Builder and
 * Source_Init will be provided.
 *
 * An additional note: non-pure virtual functions are not always required by
 * the factory builders whereas pure virtual functions are always required.  
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Interface 
{
  public:
    // usefull typedefs
    typedef std::string std_string;
    typedef std::vector<double> double_vec;
    typedef std::vector<int>    int_vec;
    typedef std::vector<std::string> string_vec;
    typedef std::vector<std::vector<int> > int_dvec;
    typedef dsxx::SP<rtt_matprops::InterpedMaterialProps> SP_Matprop;
    
  public:
    // constructor
    Interface() { /* no data to construct */ }

    // virtual constructor to make life happy down the inhertiance chain
    virtual ~Interface() { /* need a destructor for inheritance chain */ }

    // functions required by Opacity Builder

    // material data
    virtual double_vec get_density() const = 0;
    virtual double_vec get_temperature() const = 0;

    // strings that determine how to calculate the opacity and specific heats
    virtual std_string get_analytic_opacity() const = 0;
    virtual std_string get_analytic_sp_heat() const = 0;
    virtual double_vec get_kappa() const = 0;
    virtual double_vec get_kappa_thomson() const = 0;
    virtual double_vec get_specific_heat() const = 0;

    // Fleck implicitness factor
    virtual double get_implicit() const = 0;
   
    // timestep
    virtual double get_delta_t() const = 0;

    // material data
    virtual int_vec get_material_id() const = 0;
    virtual SP_Matprop get_matprops() const = 0;

    // functions required by Source_Init and Parallel_Source_Init

    // used by both parallel_source_init and source_init
    virtual double get_rad_s_tend() const = 0;
    virtual double_vec get_rad_temp() const = 0;
    virtual int get_npmax() const = 0;
    virtual int get_npnom() const = 0;
    virtual double get_dnpdt() const = 0;

    // functions defining sources that do not exist in all hosts
    virtual double_vec get_evol_ext() const { return double_vec(); }
    virtual double_vec get_rad_source() const { return double_vec(); }
    virtual string_vec get_ss_pos() const { return string_vec(); }
    virtual double_vec get_ss_temp() const { return double_vec(); }
    virtual std_string get_ss_dist() const { return std_string(); }
    virtual int get_num_global_cells() const { return int(0); }
    virtual int_vec get_global_cells() const { return int_vec(); }
    virtual int get_cycle() const { return int(0); }
    virtual double get_elapsed_t() const { return double(0); }
    virtual int get_capacity() const { return int(0); }
    virtual int_dvec get_defined_surcells() const { return int_dvec(0); }
};

} // end namespace rtt_imc

#endif                          // __imc_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Interface.hh
//---------------------------------------------------------------------------//
