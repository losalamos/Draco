//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Gray_Particle.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:30:27 2002
 * \brief  Gray_Particle class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Gray_Particle_hh__
#define __imc_Gray_Particle_hh__

#include "Particle.hh"
#include "Frequency.hh"

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Gray_Particle
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class Gray_Particle : public Particle<MT>
{
  public:

    // >>> NESTED TYPES

    /*!
     * \class Gray_Particle::Diagnostic
     * \brief Diagnostic class for tracking particle histories.
     */
    class Diagnostic : public Particle<MT>::Diagnostic
    {
      public:
	//! Constructor.
	Diagnostic(std::ostream &output_, bool detail_ = false) 
	    : Particle<MT>::Diagnostic(output_, detail_) {}
	
	// Diagnostic print functions.
	void print_xs(const Opacity<MT,Gray_Frequency> &, int) const;
    };
    
    // Friend declarations.
    friend class Diagnostic;

    // Useful typedefs.
    typedef std::vector<double>       sf_double;
    typedef rtt_rng::Sprng            Rnd_Type;
    typedef rtt_dsxx::SP<Rnd_Type>    SP_Rnd_Type;
    typedef std::string               std_string;
    typedef rtt_dsxx::SP<Diagnostic>  SP_Diagnostic;

  private:
    // >>> IMPLEMENTATION

    // Stream for implicit capture.
    inline void stream_implicit_capture(const Opacity<MT,Gray_Frequency> &,
					Tally<MT> &, double);

    // Process a collision event.
    inline void collision_event(const MT &, Tally<MT> &, double, double, 
				double);

  public:
    // Particle constructor.
    inline Gray_Particle(const sf_double &, const sf_double &, double, int,
			 Rnd_Type, double = 1, double = 1, int = BORN);

    // Unpacking constructor.
    inline Gray_Particle(const std::vector<char> &);

    // >>> TRANSPORT INTERFACE

    // IMC transport step.
    void transport(const MT &, const Opacity<MT,Gray_Frequency> &, 
		   Tally<MT> &, SP_Diagnostic = SP_Diagnostic()); 

    // >>> DIAGNOSTIC FUNCTIONS

    // Overloaded equals operator.
    bool operator==(const Gray_Particle<MT> &) const;

    //! Overloaded not equals operator.
    bool operator!=(const Gray_Particle<MT> &) const;
 
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
 * \brief Gray_Particle constructor.
 *
 * The constructor is declared inline for optimization purposes.
 */
template<class MT>
Gray_Particle<MT>::Gray_Particle(const sf_double& r_, 
				 const sf_double& omega_, 
				 double           ew_,
				 int              cell_, 
				 Rnd_Type         random_, 
				 double           frac,
				 double           tleft, 
				 int              desc)
    : Particle<MT>(r_, omega_, ew_, cell_, random_, frac, tleft, desc)
{
    // non-default particle constructor
    Check (r.size() < 4);
    Check (omega.size() == 3);
    Check (cell > 0); 
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
Gray_Particle<MT>::Gray_Particle(const std::vector<char> &packed)
{
    Require (packed.size() >= 3 * sizeof(int) + 6 * sizeof(double));

    // make an unpacker
    rtt_dsxx::Unpacker u;
    
    // set it
    u.set_buffer(packed.size(), &packed[0]);

    // unpack the spatial dimension of the particle
    int dimension = 0;
    u >> dimension;
    Check (dimension > 0 && dimension <= 3);

    // size the dimension and direction 
    r.resize(dimension);
    omega.resize(3);

    // unpack the position
    for (int i = 0; i < dimension; i++)
	u >> r[i];

    // unpack the rest of the data
    u >> omega[0] >> omega[1] >> omega[2] >> cell >> ew >> time_left
      >> fraction;
    Check (time_left >= 0.0);
    Check (fraction  >= 0.0);
    Check (cell      >  0);
    Check (ew        >= 0.0);

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
    random = new Rnd_Type(prn);
    Check (random->get_num() >= 0);
    Check (random->get_id());

    // assign the descriptor and status
    descriptor = UNPACKED;
    alive      = true;
    
    Ensure (status());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack a particle into a char stream for communication and
 * persistence. 
 */
template<class MT>
std::vector<char> Gray_Particle<MT>::pack() const
{
    Require (omega.size() == 3);

    // make a packer
    rtt_dsxx::Packer p;

    // first pack the random number state
    std::vector<char> prn = random->pack();

    // determine the size of the packed particle: 1 int for cell, + 1 int for
    // size of packed RN state + 1 int for dimension of space; dimension +
    // 6 doubles; + size of RN state chars
    int size = 3 * sizeof(int) + (r.size() + 6) * sizeof(double) + prn.size();

    // set the packed buffer
    std::vector<char> packed(size);
    p.set_buffer(size, &packed[0]);

    // pack the spatial dimension
    p << static_cast<int>(r.size());
    
    // pack the dimension
    for (int i = 0; i < r.size(); i++)
	p << r[i];
    
    // pack the rest of the data
    p << omega[0] << omega[1] << omega[2] << cell << ew << time_left
      << fraction;

    // pack the RN state
    p << static_cast<int>(prn.size());
    for (int i = 0; i < prn.size(); i++)
	p << prn[i];

    Ensure (p.get_ptr() == &packed[0] + size);
    return packed;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Stream for implicit capture.
 */
template<class MT>
void Gray_Particle<MT>::stream_implicit_capture(
    const Opacity<MT,Gray_Frequency> &xs, 
    Tally<MT>                        &tally,
    double                            distance)
{
    Check(distance >= 0.0);
    
    // exponential argument
    double argument = -xs.get_sigeffabs(cell) * distance;

    // calculate multiplicative reduction in energy-weight
    double factor = exp(argument);
    Check (factor >= minwt_frac);

    // calculate new energy weight; change in energy-weight
    double new_ew = ew * factor;
    double del_ew = ew - new_ew;

    // tally deposited energy and momentum
    tally.deposit_energy( cell, del_ew );
    tally.accumulate_momentum(cell, del_ew, omega);

    // accumulate tallies for energy-weighted path length 
    //     ewpl == energy-weighted-path-length = 
    //             int_0^d e^(-sig*x) dx =
    //             (1/sig)ew(1-e^(-sig*x)), 
    //             or if sig=0, 
    //             ewpl = ew*d.
    if (xs.get_sigeffabs(cell) > 0)
    {
	// integrate exponential from 0 to distance
	tally.accumulate_ewpl(cell, del_ew / xs.get_sigeffabs(cell) );
    }
    else if (xs.get_sigeffabs(cell) == 0)
    {
	// integrate constant
	tally.accumulate_ewpl(cell, distance * ew);
    }
    else if (xs.get_sigeffabs(cell) < 0)
    {
	Insist (0, "Effective absorption is negative!");
    }

    // update the fraction of the particle's original weight
    fraction *= factor;

    // update particle energy-weight
    ew = new_ew;

    Check(ew > 0.0);

    // Physically transport the particle
    stream( distance ); 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a collision.
 */
template<class MT> 
void Gray_Particle<MT>::collision_event(
    const MT  &mesh, 
    Tally<MT> &tally, 
    double     prob_scatter, 
    double     prob_thomson_scatter, 
    double     prob_abs)
{
    Check (mesh.num_cells() == tally.num_cells());

    // get a random number
    double rand_selector = random->ran();
    
    if (rand_selector < prob_scatter) 
    { 
	// Scatter
	if (rand_selector < prob_thomson_scatter)
	{ 
	    // Thomson scatter
	    descriptor = THOM_SCATTER;
	    tally.accum_n_thomscat();
	}
	else
	{ 
	    // Effective scatter
	    descriptor = EFF_SCATTER;
	    tally.accum_n_effscat();
	}
	
	// accumulate momentum from before the scatter
	tally.accumulate_momentum(cell, ew, omega);
	
	// scatter the particle -- update direction cosines (we really should
	// call effective_scatter for effective scatters because the particle
	// is re-emitted so the direction should be sampled from rest, since
	// everything is isotropic we will wait to do this later (so as not
	// to hose all our regression tests); it has no effect on the
	// practical outcome of the transport
	scatter( mesh );
	
	// accumulate momentum from after the scatter
	tally.accumulate_momentum(cell, -ew, omega);
	
    }
    else if (rand_selector < prob_scatter + prob_abs)
    { 
	// Absorption
	
	// tally absorption data
	tally.deposit_energy( cell, ew );
	tally.accum_n_killed();
	tally.accum_ew_killed( ew );
	tally.accumulate_momentum(cell, ew, omega);

	// set the descriptor and particle status
	descriptor = KILLED; 
	alive      = false;
    }
    else
    {
	Insist(0,"D'oh! Transport could not pick a random event!");
    }   
}

} // end namespace rtt_imc

#endif                          // __imc_Gray_Particle_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Gray_Particle.hh
//---------------------------------------------------------------------------//
