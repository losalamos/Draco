//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Source_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Dec  8 14:35:42 1999
 * \brief  Source_Builder header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Source_Builder_hh__
#define __imc_Source_Builder_hh__

#include "Particle.hh"
#include "Global.hh"
#include "Source.hh"
#include "Mat_State.hh"
#include "Opacity.hh"
#include "mc/Topology.hh"
#include "mc/Particle_Stack.hh"
#include "mc/Parallel_Data_Operator.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>
#include <cmath>
#include <typeinfo>

namespace rtt_rng
{

// Forward declarations.
class Rnd_Control;
class Sprng;

} // end of namespace rtt_rng

namespace rtt_mc
{

// Forward declarations.
class Comm_Patterns;

}

namespace rtt_imc
{

class Gray_Frequency;
class Multigroup_Frequency;
template<class MT> class Gray_Particle;
template<class MT> class Multigroup_Particle;
 
//===========================================================================//
/*!
 * \class Source_Builder
 *
 * \brief Source_Builder abstract base class for building IMC sources.
 *
 * Source_Builder is an abstract base class used to create source builder
 * inheritance hierarchies.  In short, Source_Builder builds an IMC source on
 * processor.  The derived classes of Source_Builder are: 
 *
 * \arg rtt_imc::Rep_Source_Builder for full replication sources; \arg
 * rtt_imc::DD_Source_Builder for full domain decomposition sources; \arg
 * rtt_imc::Serial_Source_Builder for serial sources.
 *
 * Rep_Source_Builder defaults to serial when run on one processor.
 * Serial_Source_Builder is simply an optimized version of Rep_Source_Builder
 * for 1 processor (non-communication) compilations.  Thus, we have full
 * replication, full DD, DD/replication, and serial derived classes.
 * Functionality common to all types of source builders exists in the base
 * class as non-virtual functions.
 *
 * Source_Builder and Source_Builder hierarchies are templated on mesh and
 * particle type.  By default, the particle type is rtt_imc::Particle.
 *
 * Source_Builder classes use the rtt_mc::Topology classes to determine where
 * things are and where they need to go.
 *
 * Source_Builder is used to build smart pointers to rtt_imc::Source on each
 * processor.  Additionally, the Source_Builder class, in the process of
 * making sources, calculates global source information that can be used for
 * edits.  Access functions are provided for this data as regular and virtual
 * functions that are defined in the derived classes 
 *
 * All protected data in the Source_Builder base class is local to a
 * processor.  Global quantities are calculated by using the
 * Parallel_Data_Operations::local_to_global mapper and specifying whether
 * the data to be mapped is Data_Replicated or Data_Decomposed.
 * Specifically, when every replicated cell has exactly the same value, the
 * data is replicated.  When replicated cells have a portion of the global
 * value, the data is said to be decomposed.  If no cells are replicated,
 * i.e., the topology parallel scheme is full domain decomposition, all data
 * is distributed (Data_Distibuted).  In the case of full replication,
 * deterministic quantities such as volume emission energy and surface source
 * energy are data-replicated, and stochastic quantities such as the
 * surviving census particle energies are data-decomposed.
 *
 * The following relationships apply to topology and data decomposition for
 * locally calculated data in "replication" topologies: \arg local data +
 * Data_Replication = global data \arg local data + Data_Decomposed = local
 * data
 *
 * Data_Distributed does not apply to "replication".  In full "DD" all data
 * is Data_Distributed, thus the local data on each processor is inherently
 * global for that cell.  A caveat is that global data fields in "DD"
 * topologies may not contain all the global data elements on a given
 * processor.  It may only contain the data that lives on that processor,
 * ie. the data that is mapped to local cells.
 *
 * Random number stream ID fields are set using the rtt_rng::rn_stream
 * variable.  This variable is \b defined to be global in \b all topologies.
 * For some topologies ("DD") this requires communication after random number
 * IDs are assigned to set and to insure that the rtt_rng::rn_stream value is
 * the same on all processors. 
 *
 * Accessors to data that live in the Source_Builder base class returns the
 * values that have been calculated on processor.  It is up to the client to
 * do any topology-dependent mechanics required to get global (problem-wide)
 * data.  The data could be global or local depending upon the topology.
 * Virtual functions are used to access source data that live in derived
 * classes.  Accordingly, this data is also global or local depending upon
 * the topology.  
 */
/*!
 * \example imc/test/tstSource_Builder.cc
 *
 * Example usage of the Source_Builder class.  The derived classes
 * Rep_Source_Builder and DD_Source_Builder are tested.
 */
// revision history:
// -----------------
// 0) original
// 1) 19 Jun 2000 : added function reset_ew_in_census to allow for post-comb
//                  census energy-weight readjustment.
// 2) 03-JUL-00   : made get_initial_census_energy a virtual function because
//                  in DD topologies ecentot is a local value and in full
//                  replication it is a global value
// 3) 24 Aug 2000 : added capability for the random number stream ID's to
//                  wrap around 2e9.
// 4) 31-JUL-2001 : changed mod_with_2e9 to INTEGER_MODULO_1E9 from rtt_mc
// 4) 07-JAN-2002 : moved constructor to header file so that automatic
//                  instantiation will work; updated to work with new
//                  Particle_Stack 
//===========================================================================//

template<class MT, class FT, class PT>
class Source_Builder 
{
  public:
    // typedefs used in the inheritance chain
    typedef rtt_dsxx::SP<Source<MT,FT,PT> >       SP_Source;
    typedef rtt_mc::Particle_Containers<PT>       Containers;
    typedef typename Containers::Census           Census;
    typedef rtt_dsxx::SP<Census>                  SP_Census;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >         SP_Opacity;
    typedef rtt_dsxx::SP<Mat_State<MT> >          SP_Mat_State;
    typedef rtt_dsxx::SP<MT>                      SP_Mesh;
    typedef rtt_dsxx::SP<rtt_rng::Rnd_Control>    SP_Rnd_Control;
    typedef rtt_dsxx::SP<rtt_mc::Topology>        SP_Topology;
    typedef rtt_dsxx::SP<rtt_mc::Comm_Patterns>   SP_Comm_Patterns;
    typedef std::vector<int>                      sf_int;
    typedef std::vector<std::vector<int> >        vf_int;
    typedef std::vector<double>                   sf_double;
    typedef std::vector<std::string>              sf_string;
    typedef std::string                           std_string;
    typedef typename MT::template CCSF<double>    ccsf_double;
    typedef typename MT::template CCSF<int>       ccsf_int;
    typedef typename MT::template CCVF<double>    ccvf_double;

  private:
    // >>> BASE CLASS DATA

    // Interface data.

    //! Elapsed time that problem has been run.
    double elapsed_t;

    //! Vector of local external volume source energies.
    sf_double evol_ext;

    //! External radiation source duration.
    double rad_s_tend;

    //! External radiation source energies.
    sf_double rad_source;

    //! Initial radiation temperature, which determine initial census.
    sf_double rad_temp;

    //! Position (and number) of surface sources.
    sf_string ss_pos;

    //! Temperatures of surface sources.
    sf_double ss_temp;

    //! List of local cells that have surface sources.
    vf_int defined_surcells;

    //! Descriptor for each surface source.
    sf_string ss_desc;

    // >>> BASE CLASS IMPLEMENTATION MEMBER FUNCTIONS

    // Calculate volume emission source energies.
    void calc_evol(const Mat_State<MT> &, const Opacity<MT,FT> &);

    // Calculate surface source energies.
    void calc_ess(const Opacity<MT,FT> &);

  private:
    // >>> SPECIALIZATIONS ON FREQUENCY
    
    // Gray_Frequency partial specialization of calc_evol_net.
    inline void calc_evol_net(
	rtt_imc::global::Type_Switch<Gray_Frequency>,
	const Mat_State<MT> &,
	const Opacity<MT,Gray_Frequency> &);
    
    // Multigroup_Frequency partial specialization of calc_evol_net.
    inline void calc_evol_net(
	rtt_imc::global::Type_Switch<Multigroup_Frequency>,
	const Mat_State<MT> &, 
	const Opacity<MT,Multigroup_Frequency> &);
    
    // Calculate probability of straight Planckian emission for
    // Gray_Frequency.
    template<class Stop_Explicit_Instantiation>
    inline void calc_prob_Planck_emission(
	rtt_imc::global::Type_Switch<Gray_Frequency>, double, int);

    // Calculate probability of straight Planckian emission for
    // Multigroup_Frequency.
    template<class Stop_Explicit_Instantiation>
    inline void calc_prob_Planck_emission(
	rtt_imc::global::Type_Switch<Multigroup_Frequency>, double, int);

    // Calculate particle for Gray_Particle type.
    template<class Stop_Explicit_Instantiation>
    rtt_dsxx::SP<Gray_Particle<MT> > make_particle(
	rtt_imc::global::Type_Switch<Gray_Frequency>,
	rtt_imc::global::Type_Switch<Gray_Particle<MT> >,
	Gray_Frequency, const double, const sf_double &,
	const sf_double &, const double, int, 
	const rtt_rng::Sprng &);

    // Calculate particle for Multigroup_Particle type.
    template<class Stop_Explicit_Instantiation>
    rtt_dsxx::SP<Multigroup_Particle<MT> > make_particle(
	rtt_imc::global::Type_Switch<Multigroup_Frequency>,
	rtt_imc::global::Type_Switch<Multigroup_Particle<MT> >,
	Multigroup_Frequency, const double, const sf_double &,
	const sf_double &, const double, int, 
	const rtt_rng::Sprng &);

  protected:
    // >>> DATA USED BY ALL SOURCE BUILDERS

    // Interface data.

    //! Requested number of source particles.
    int npnom;

    //! Maximum number of source particles.
    int npmax;

    //! Rate of change of number of source particles per shake.
    double dnpdt;

    //! Problem cycle.
    int cycle;

    //! Timestep.
    double delta_t;

    //! Surface source angular distribution.
    std_string ss_dist;
    
    //! Local (on-processor) census.
    SP_Census census;

    // Data through constructor arguments.

    //! Problem topology.
    SP_Topology topology;

    //! Parallel_Data_Operator for parallel data operations.
    rtt_mc::Parallel_Data_Operator parallel_data_op;

    // Source class data.

    //! Number of particles that we want for this cycle.
    int npwant;

    //! Local field of census energies.
    ccsf_double ecen;

    //! Local (initial) census total energy.
    double ecentot;

    //! Local field of census energy weights.
    ccsf_double ew_cen;

    //! Local field of volume emission energies.
    ccsf_double evol;

    //! Local field of net volume emission energies.
    ccsf_double evol_net;
    
    //! Local field of external material volume source energies.
    ccsf_double mat_vol_src;

    //! Local volume emission total energy.
    double evoltot;

    //! Local field of volume emission energy weights.
    ccsf_double ew_vol;

    //! Local external material volume source total energy.
    double mat_vol_srctot;

    //! Local field of surface source energies.
    ccsf_double ess;

    //! Local field of surface source faces.
    ccsf_int ss_face_in_cell;

    //! Local surface source total energy.
    double esstot;

    //! Local field of surface source energy weights.
    ccsf_double ew_ss;

    //! Local starting random number stream ID field for volume emission.
    ccsf_int volrn;

    //! Local starting random number stream ID field for surface source.
    ccsf_int ssrn;

    //! Frequency-dependent sampling of evol_net.
    Frequency_Sampling_Data<MT,FT> freq_samp_data;

    // >>> IMPLEMENTATION INHERITANCE

    // Calculate source energies for volume emission and surface source.
    void calc_source_energies(const Mat_State<MT> &, const Opacity<MT,FT> &);

    // Calculate the initial census energy.
    void calc_initial_ecen(const Opacity<MT,FT> &);

    // Calculate the number of source particles.
    void calc_num_src_particles(const double, const ccsf_double &, 
				ccsf_int &, int &); 

    // Write the local census on each processor.
    void write_initial_census(SP_Mesh, SP_Rnd_Control, const FT &,
			      const Mat_State<MT> &, const ccsf_int &, 
			      const int &, const ccsf_int &);

    // Comb local census.
    void comb_census(SP_Rnd_Control, ccsf_int &, int &, double &, SP_Census,
		     ccsf_int &);

    // Reset the energy-weight of each particle in the local census.
    void reset_ew_in_census(const int &, double &, const double &);

  public:
    // Constructor.
    template<class IT>
    Source_Builder(rtt_dsxx::SP<IT>, SP_Mesh, SP_Topology);

    // virtual destructor for correct behavior in inheritance chain
    virtual ~Source_Builder() {/*...*/}

    //! Build rtt_dsxx::SP to source.
    virtual SP_Source build_Source(SP_Mesh, SP_Mat_State, SP_Opacity,
				   SP_Rnd_Control, SP_Comm_Patterns) = 0;

    //! Build an initial census from initial radiation temperature.
    virtual void calc_initial_census(SP_Mesh, SP_Mat_State, SP_Opacity,
				     SP_Rnd_Control) = 0;

    //! Get the number of cells accessed by the Source_Builder on processor.
    int num_cells() const { return topology->num_cells(C4::node()); }

    // >>> ACCESSOR FOR CENSUS

    //! Get a rtt_dsxx::SP to the census.
    SP_Census get_census() const { return census; }

    // >>> ACCESSOR SERVICES TO GET SOURCE INITIALIZATION DATA

    // These functions are primarily for edits and the problem variables on
    // processor.  The client will be required to use the data in an
    // appropriate fashion based on the problem topology, ie. data could be
    // either global or local depending upon the topology.

    //! Get total volume emission energy on processor.
    double get_evoltot() const { return evoltot; } 

    //! Get global material volume emission source energy/cell on processor.
    double get_mat_vol_src(int cell) const { return mat_vol_src(cell); } 

    //! Get total material volume emission source energy on processor.
    double get_mat_vol_srctot() const { return mat_vol_srctot; }

    //! Get net volume emission energy per cell on processor.
    double get_evol_net(int cell) const { return evol_net(cell); }

    //! Get total surface source energy on processor.
    double get_esstot() const { return esstot; }

    // Topology dependent data.  This data could be either global or local
    // depending upon the topology.  The data that is returned lives in the
    // derived class.

    //! Get total initial census energy on processor - topology dependent.
    virtual double get_initial_census_energy() const = 0;

    //! Get energy loss in volume emission - topology dependent.
    virtual double get_eloss_vol() const = 0;

    //! Get total number of volume emission particles - toplogy dependent.
    virtual int get_nvoltot() const = 0;

    //! Get number of volume emission particles / cell - topology dependent.
    virtual int get_nvol(int) const = 0;
    
    //! Get energy loss in surface source - topology dependent.
    virtual double get_eloss_ss() const = 0;

    //! Get total number of surface source particles - topology dependent.
    virtual int get_nsstot() const = 0;

    //! Get energy loss in census - topology dependent.
    virtual double get_eloss_cen() const = 0;
    
    //! Get total number of post-comb census particles - topology dependent.
    virtual int get_ncentot() const = 0;
};

//---------------------------------------------------------------------------//
// TEMPLATE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Source_Builder base class constructor.
 *
 * The Source_Builder base class constructor gets data from a run-time
 * interface and a valid Topology to construct its data.  The data in
 * Source_Builder is divided into two basic types: \arg interface data comes
 * from the interface and is used to calculate source fields, \arg source
 * data fields are the products of the Source_Builder process.  The interface
 * data must be given to Source_Builder locally.  That is, the data should be
 * dimensioned to the local (on-processor) mesh size.  All fields in
 * Source_Builder are dimensioned to the local mesh.
 *
 * When the constructor is called the following requirements must be met:
 * \arg a mesh must exist and be properly sized, \arg the rtt_rng::rn_stream
 * must be the same on all processors.
 *
 * \param interface rtt_dsxx::SP to a valid run-time interface
 * \param mesh rtt_dsxx::SP to a local mesh object
 * \param top rtt_dsxx::SP to a topology object 
 */
template<class MT, class FT, class PT>
template<class IT>
Source_Builder<MT,FT,PT>::Source_Builder(rtt_dsxx::SP<IT> interface,
					 SP_Mesh          mesh, 
					 SP_Topology      top)
    : elapsed_t(interface->get_elapsed_t()), 
      evol_ext(interface->get_evol_ext()),
      rad_s_tend(interface->get_rad_s_tend()),
      rad_source(interface->get_rad_source()),
      rad_temp(interface->get_rad_temp()),
      ss_pos(interface->get_ss_pos()),
      ss_temp(interface->get_ss_temp()),
      defined_surcells(interface->get_defined_surcells()),
      ss_desc(interface->get_ss_desc()),
      npnom(interface->get_npnom()),
      npmax(interface->get_npmax()),
      dnpdt(interface->get_dnpdt()),
      cycle(interface->get_cycle()), 
      delta_t(interface->get_delta_t()), 
      ss_dist(interface->get_ss_dist()),
      topology(top), 
      parallel_data_op(top),
      census(interface->get_census()),
      npwant(0),
      ecen(mesh),
      ew_cen(mesh),
      ecentot(0),
      evol(mesh),
      ew_vol(mesh),
      evol_net(mesh),
      mat_vol_src(mesh),
      evoltot(0),
      mat_vol_srctot(0),
      ess(mesh),
      ew_ss(mesh),
      ss_face_in_cell(mesh),
      esstot(0),
      volrn(mesh),
      ssrn(mesh),
      freq_samp_data(mesh)
{
    using rtt_mc::global::min;

    Require(mesh);
    Require(mesh->num_cells() == topology->num_cells(C4::node()));
    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));

    // modulo the rn_stream with 1e9 so that, when we get to more than
    // 1e9 particles (each with its own rn_stream, numbered 0 to 1e9-1), the
    // rn_stream starts back with rn_stream=0.  The rng package is still
    // limited to rnstream < numgen, so the 1e9 wrap-around is moot unless
    // numgen is 1e9 or higher (int size limiting + spawn change).
    rtt_rng::rn_stream = INTEGER_MODULO_1E9(rtt_rng::rn_stream);

    int num_cells = mesh->num_cells();

    // calculate the desired number of source particles
    npwant = min(npmax, static_cast<int>(npnom + dnpdt * elapsed_t)); 

    Ensure(evol_ext.size() == num_cells);
    Ensure(rad_source.size() == num_cells);
    Ensure(rad_temp.size() == num_cells || rad_temp.size() == 0);
    Ensure(ss_pos.size() == ss_temp.size());
    Ensure(ss_pos.size() == defined_surcells.size());
    Ensure(npwant > 0);
}

//---------------------------------------------------------------------------//
// PARTIAL SPECIALIZATIONS BASED ON FREQUENCY (PRIVATE)
//---------------------------------------------------------------------------//
// Gray_Frequency partial specialization of calc_evol_net.

template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::calc_evol_net(
    rtt_imc::global::Type_Switch<Gray_Frequency>,
    const Mat_State<MT>              &state,
    const Opacity<MT,Gray_Frequency> &opacity)
{
    using rtt_imc::global::Type_Switch;
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using std::pow;

    Check (typeid(FT) == typeid(Gray_Frequency));

    for (int cell = 1; cell <= evol.size(); cell++)
	evol_net(cell) = opacity.get_fleck(cell) *
	    opacity.get_sigma_abs(cell) * a * c * pow(state.get_T(cell), 4) *
	    evol.get_Mesh().volume(cell) * delta_t;
}

//---------------------------------------------------------------------------//
// Multigroup_Frequency partial specialization of calc_evol_net.

template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::calc_evol_net(
    rtt_imc::global::Type_Switch<Multigroup_Frequency>,
    const Mat_State<MT>                    &state,
    const Opacity<MT,Multigroup_Frequency> &opacity)
{
    using rtt_imc::global::Type_Switch;
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using std::pow;

    Check (typeid(FT) == typeid(Multigroup_Frequency));

    for (int cell = 1; cell <= evol.size(); cell++)
	evol_net(cell) = opacity.get_fleck(cell) *
	    opacity.get_integrated_sigma_times_Planck(cell) * a * c *
	    pow(state.get_T(cell), 4) * evol.get_Mesh().volume(cell) *
	    delta_t;
}

//---------------------------------------------------------------------------//
// Calculate probability of straight Planckian emission for Gray_Frequency

template<class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
void Source_Builder<MT,FT,PT>::calc_prob_Planck_emission(
    rtt_imc::global::Type_Switch<Gray_Frequency>,
    double evol_add,
    int    cell)
{
    // NULL OP FOR GRAY
}

//---------------------------------------------------------------------------//
// Calculate probability of straight Planckian emission for
// Multigroup_Frequency.

template<class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
void Source_Builder<MT,FT,PT>::calc_prob_Planck_emission(
    rtt_imc::global::Type_Switch<Multigroup_Frequency>,
    double evol_add,
    int    cell)
{
    using rtt_imc::global::Type_Switch;

    Check (typeid(FT) == typeid(Multigroup_Frequency));

    // calculate the probability of straight Planckian emission in a cell
    freq_samp_data.prob_of_straight_Planck_emission(cell) = evol_add / 
	evol(cell);
}

} // end namespace rtt_imc

#endif                          // __imc_Source_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source_Builder.hh
//---------------------------------------------------------------------------//
