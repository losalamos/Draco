//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Multigroup_Particle.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:31:05 2002
 * \brief  Multigroup_Particle class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_imc_Multigroup_Particle_HH
#define RTT_imc_Multigroup_Particle_HH

#include "Particle.hh"
#include "Random_Walk.hh"
#include "Extrinsic_Surface_Tracker.hh"
#include "mc/Sampler.hh"

namespace rtt_imc
{
 
// Forward declarations.
class Multigroup_Frequency;

//===========================================================================//
/*!
 * \class Multigroup_Particle
 *
 * \brief Particle for doing multigroup IMC transport.
 */
// revision history:
// -----------------
// 0) original
// 1) 10-FEB-03 : removed #define's for Base class scoping; added real
//                scoping
// 2) 19-MAY-03 : added Random_Walk to transport
// 3) 29-JUL-03 : implemented Random_Walk in transport
// 
//===========================================================================//

template<class MT>
class Multigroup_Particle : public Particle<MT>
{
  public:

    // >>> NESTED TYPES

    /*!
     * \brief Diagnostic class for tracking multigroup particle histories.
     */
    class Diagnostic : public Particle<MT>::Diagnostic
    {
      public:
	//! Constructor.
	Diagnostic(std::ostream &output_, bool detail_ = false) 
	    : Particle<MT>::Diagnostic(output_, detail_) {}
	
	// Diagnostic print functions.
	void print_xs(const Opacity<MT,Multigroup_Frequency> &, int, int) const;
    };
    
    // Friend declarations.
    friend class Diagnostic;

    // Useful typedefs.
    typedef std::vector<double>                     sf_double;
    typedef rtt_rng::Sprng                          Rnd_Type;
    typedef rtt_dsxx::SP<Rnd_Type>                  SP_Rnd_Type;
    typedef std::string                             std_string;
    typedef rtt_dsxx::SP<Diagnostic>                SP_Diagnostic;
    typedef Opacity<MT,Multigroup_Frequency>        MG_Opacity;
    typedef rtt_dsxx::SP<Random_Walk<MT> >          SP_Random_Walk;
    typedef rtt_dsxx::SP<Extrinsic_Surface_Tracker> SP_Surface_tracker;

  private:
    // Typedef for base class scoping.
    typedef Particle<MT> Base;

  private:
    // >>> DATA

    // Group index [1,N_groups] where 1 is the lowest-frequency group.
    int group_index;

  private:
    // >>> IMPLEMENTATION

    // Transport a particle without hybrid methods.
    void straight_transport(const MT &, const MG_Opacity &, Tally<MT> &, 
			    SP_Surface_tracker, SP_Diagnostic);

    // Transport a particle with random walk.
    void rw_transport(const MT &, const MG_Opacity &, Tally<MT> &, 
		      SP_Random_Walk, SP_Surface_tracker, SP_Diagnostic);

    // Process a collision event.
    void collision_event(const MT &, Tally<MT> &, const MG_Opacity &,
			 double, double, double);

    // Do an effective scatter
    inline void effective_scatter(const MT &, const MG_Opacity &);

  public:
    // Particle constructor.
    inline Multigroup_Particle(const sf_double &, const sf_double &, double,
			       int, Rnd_Type, int, double = 1, double = 1, 
			       int = Base::BORN);

    // Unpacking constructor.
    Multigroup_Particle(const std::vector<char> &);

    // >>> TRANSPORT INTERFACE

    // IMC transport step.
    void transport(const MT &, const MG_Opacity &, Tally<MT> &, 
		   SP_Random_Walk = SP_Random_Walk(),
		   SP_Surface_tracker = SP_Surface_tracker(),
		   SP_Diagnostic = SP_Diagnostic()); 

    // >>> ACCESSORS

    //! Get the particle frequency group index.
    int get_group_index() const { return group_index; }

    // >>> DIAGNOSTIC FUNCTIONS

    // Print the particle's state.
    void print(std::ostream &) const;

    // Overloaded equals operator.
    bool operator==(const Multigroup_Particle<MT> &) const;

    //! Overloaded not equals operator.
    bool operator!=(const Multigroup_Particle<MT> &) const;
 
    // >>> PACKING FUNCTIONALITY
    
    // Pack function
    std::vector<char> pack() const;

    // Get the size of the packed particle.
    static int get_packed_particle_size(int, const rtt_rng::Rnd_Control &);
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Multigroup_Particle constructor.
 *
 * The constructor is declared inline for optimization purposes.
 */
template<class MT>
Multigroup_Particle<MT>::Multigroup_Particle(const sf_double& r_, 
					     const sf_double& omega_, 
					     double           ew_,
					     int              cell_, 
					     Rnd_Type         random_, 
					     int              group_index_,
					     double           frac,
					     double           tleft, 
					     int              desc)
    : Particle<MT>(r_, omega_, ew_, cell_, random_, frac, tleft, desc),
      group_index(group_index_)
{
    // non-default particle constructor
    Check (Base::r.size() < 4);
    Check (Base::omega.size() == 3);
    Check (Base::cell > 0); 
    Check (group_index > 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do an effective scatter.
 *
 * Sample a new frequency, for mulitgroup particles in the Fleck and
 * Cummings method represents an absorption followed by emission.  Thus,
 * we need to sample a new frequency group on the re-emission.  Also, we need
 * to sample a new particle direction.
 */
template<class MT>
void Multigroup_Particle<MT>::effective_scatter(const MT         &mesh,
						const MG_Opacity &opacity)
{
    using rtt_mc::sampler::sample_bin_from_discrete_cdf;
    
    // sample new group index (add 1 since return value is in [0,G-1])
    group_index = sample_bin_from_discrete_cdf(
	*Base::random, opacity.get_emission_group_cdf(Base::cell)) + 1;
    Check (group_index > 0);
    Check (group_index <= opacity.get_Frequency()->get_num_groups());
    
    // sample particle direction for re-emitted particle
    Base::omega = mesh.get_Coord().sample_dir("isotropic", *Base::random);
}

} // end namespace rtt_imc

#endif                          // RTT_imc_Multigroup_Particle_HH

//---------------------------------------------------------------------------//
//                              end of imc/Multigroup_Particle.hh
//---------------------------------------------------------------------------//
