//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Multigroup_Particle.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:31:05 2002
 * \brief  Multigroup_Particle class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_imc_Multigroup_Particle_HH
#define RTT_imc_Multigroup_Particle_HH

#include "Particle.hh"
#include "Random_Walk.hh"
#include "mc/Sampler.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cmath>

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
// 
//===========================================================================//

template<class MT>
class Multigroup_Particle : public Particle<MT>
{
  public:
    // Useful typedefs.
    typedef std::vector<double>              sf_double;
    typedef rtt_rng::Sprng                   Rnd_Type;
    typedef rtt_dsxx::SP<Rnd_Type>           SP_Rnd_Type;
    typedef std::string                      std_string;
    typedef Opacity<MT,Multigroup_Frequency> MG_Opacity;
    typedef rtt_dsxx::SP<Random_Walk<MT> >   SP_Random_Walk;

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
	void print_xs(const MG_Opacity &, int, int) const;
    };

    // Diagnostic typedef.
    typedef rtt_dsxx::SP<Diagnostic> SP_Diagnostic;
    
    // Friend declarations.
    friend class Diagnostic;

  private:
    // Typedef for base class scoping.
    typedef Particle<MT> Base;

  private:
    // >>> DATA

    // Group index [1,N_groups] where 1 is the lowest-frequency group.
    int group_index;

  private:
    // >>> IMPLEMENTATION

    // Process a collision event.
    inline void collision_event(const MT &, Tally<MT> &, const MG_Opacity &,
				double, double, double);

    // Do an effective scatter
    inline void effective_scatter(const MT &, const MG_Opacity &);

  public:
    // Particle constructor.
    inline Multigroup_Particle(const sf_double &, const sf_double &, double,
			       int, Rnd_Type, int, double = 1, double = 1, 
			       int = Base::BORN);

    // Unpacking constructor.
    inline Multigroup_Particle(const std::vector<char> &);

    // >>> TRANSPORT INTERFACE

    // IMC transport step.
    void transport(const MT &, const MG_Opacity &, 
		   Tally<MT> &, SP_Random_Walk = SP_Random_Walk(),
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
    inline std::vector<char> pack() const;

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
 * \brief Unpacking constructor.
 *
 * This constructor is used to unpack a particle that has been packed with
 * the Particle::pack() function.
 *
 * It is declared inline for optimization.
 */
template<class MT>
Multigroup_Particle<MT>::Multigroup_Particle(const std::vector<char> &packed)
{
    Require (packed.size() >= 4 * sizeof(int) + 6 * sizeof(double));

    // make an unpacker
    rtt_dsxx::Unpacker u;
    
    // set it
    u.set_buffer(packed.size(), &packed[0]);

    // unpack the spatial dimension of the particle
    int dimension = 0;
    u >> dimension;
    Check (dimension > 0 && dimension <= 3);

    // size the dimension and direction 
    Base::r.resize(dimension);
    Base::omega.resize(3);

    // unpack the position
    for (int i = 0; i < dimension; i++)
	u >> Base::r[i];

    // unpack the rest of the data
    u >> Base::omega[0] >> Base::omega[1] >> Base::omega[2] >> Base::cell 
      >> Base::ew >> Base::time_left >> Base::fraction >> group_index;
    Check (Base::time_left   >= 0.0);
    Check (Base::fraction    >= 0.0);
    Check (Base::cell        >  0);
    Check (Base::ew          >= 0.0);
    Check (group_index >  0)

    // get the size of the RN state
    int size_rn = 0;
    u >> size_rn;
    Check (size_rn > 0);

    // make a packed rn vector
    std::vector<char> prn(size_rn);

    // unpack the rn state
    for (int i = 0; i < size_rn; i++)
	u >> prn[i];

    // rebuild the rn state
    Base::random = new Rnd_Type(prn);
    Check (Base::random->get_num() >= 0);
    Check (Base::random->get_id());

    // assign the descriptor and status
    Base::descriptor = Base::UNPACKED;
    Base::alive      = true;
    
    Ensure (Base::status());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack a particle into a char stream for communication and
 * persistence. 
 */
template<class MT>
std::vector<char> Multigroup_Particle<MT>::pack() const
{
    Require (Base::omega.size() == 3);

    // make a packer
    rtt_dsxx::Packer p;

    // first pack the random number state
    std::vector<char> prn = Base::random->pack();

    // determine the size of the packed particle: 1 int for cell, + 1 int for
    // size of packed RN state + 1 int for dimension of space + 1 int for
    // group index; dimension + 6 doubles; + size of RN state chars
    int size = 4 * sizeof(int) + (Base::r.size() + 6) * sizeof(double) +
	prn.size();

    // set the packed buffer
    std::vector<char> packed(size);
    p.set_buffer(size, &packed[0]);

    // pack the spatial dimension
    p << static_cast<int>(Base::r.size());
    
    // pack the dimension
    for (int i = 0; i < Base::r.size(); i++)
	p << Base::r[i];
    
    // pack the rest of the data
    p << Base::omega[0] << Base::omega[1] << Base::omega[2] << Base::cell 
      << Base::ew << Base::time_left << Base::fraction << group_index;

    // pack the RN state
    p << static_cast<int>(prn.size());
    for (int i = 0; i < prn.size(); i++)
	p << prn[i];

    Ensure (p.get_ptr() == &packed[0] + size);
    return packed;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a collision.
 */
template<class MT> 
void Multigroup_Particle<MT>::collision_event(
    const MT         &mesh, 
    Tally<MT>        &tally, 
    const MG_Opacity &opacity,
    double            prob_scatter, 
    double            prob_thomson_scatter, 
    double            prob_abs)
{
    Check (mesh.num_cells() == tally.num_cells());

    // get a random number
    double rand_selector = Base::random->ran();
    
    if (rand_selector < prob_scatter) 
    { 
	// accumulate momentum from before the scatter
	tally.accumulate_momentum(Base::cell, Base::ew, Base::omega);
	
	// Scatter
	if (rand_selector < prob_thomson_scatter)
	{ 
	    // Thomson scatter
	    Base::descriptor = Base::THOM_SCATTER;
	    tally.accum_n_thomscat();
	
	    // scatter the particle -- update direction cosines
	    Base::scatter(mesh);
	}
	else
	{ 
	    // Effective scatter
	    Base::descriptor = Base::EFF_SCATTER;
	    tally.accum_n_effscat();

	    // scatter the particle -- update direction cosines
	    effective_scatter(mesh, opacity);
	}
	
	// accumulate momentum from after the scatter
	tally.accumulate_momentum(Base::cell, -Base::ew, Base::omega);
    }
    else if (rand_selector < prob_scatter + prob_abs)
    { 
	// Absorption
	
	// tally absorption data
	tally.deposit_energy(Base::cell, Base::ew);
	tally.accum_n_killed();
	tally.accum_ew_killed(Base::ew);
	tally.accumulate_momentum(Base::cell, Base::ew, Base::omega);

	// set the descriptor and particle status
	Base::descriptor = Base::KILLED; 
	Base::alive      = false;
    }
    else
    {
	Insist(0,"D'oh! Transport could not pick a random event!");
    }   
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
