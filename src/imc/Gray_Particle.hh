//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Gray_Particle.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 29 13:30:27 2002
 * \brief  Gray_Particle class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_imc_Gray_Particle_HH
#define RTT_imc_Gray_Particle_HH

#include "Particle.hh"
#include "Random_Walk.hh"
#include "Extrinsic_Surface_Tracker.hh"

namespace rtt_imc
{

// Forward declarations
class Gray_Frequency;

//===========================================================================//
/*!
 * \class Gray_Particle
 *
 * \brief Particle for doing gray IMC transport.
 */
// revision history:
// -----------------
// 0) original
// 1) 10-FEB-03 : removed #define's for Base class scoping; added real
//                scoping
// 2) 18-FEB-03 : added random walk version of transport function
// 3) 24-JUN-03 : updated for new sub tally concept in Tally<MT>
// 
//===========================================================================//

template<class MT>
class Gray_Particle : public Particle<MT>
{
  public:

    // >>> NESTED TYPES

    /*!
     * \brief Diagnostic class for tracking gray particle histories.
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
    typedef std::vector<double>            sf_double;
    typedef rtt_rng::Sprng                 Rnd_Type;
    typedef rtt_dsxx::SP<Rnd_Type>         SP_Rnd_Type;
    typedef std::string                    std_string;
    typedef rtt_dsxx::SP<Diagnostic>       SP_Diagnostic;
    typedef rtt_dsxx::SP<Random_Walk<MT> > SP_Random_Walk;
    typedef rtt_dsxx::SP<Extrinsic_Surface_Tracker>  SP_Surface_tracker;

  private:
    // Typedef for base class scoping.
    typedef Particle<MT> Base;

  private:
    // >>> IMPLEMENTATION

    // Transport a particle without hybrid methods.
    void straight_transport(const MT &, const Opacity<MT,Gray_Frequency> &,
			    Tally<MT> &, SP_Surface_tracker, SP_Diagnostic);

    // Transport a particle with random walk.
    void rw_transport(const MT &, const Opacity<MT,Gray_Frequency> &,
		      Tally<MT> &, SP_Random_Walk, SP_Surface_tracker,
		      SP_Diagnostic);

    // Process a collision event.
    void collision_event(const MT &, Tally<MT> &, double, double, double);

    // Process a random walk particle.
    void random_walk_event(double, Tally<MT> &,
			   const Opacity<MT,Gray_Frequency> &);

  public:
    // Particle constructor.
    inline Gray_Particle(const sf_double &, const sf_double &, double, int,
			 Rnd_Type, double = 1, double = 1, int = Base::BORN);

    // Unpacking constructor.
    Gray_Particle(const std::vector<char> &);

    // >>> TRANSPORT INTERFACE 

    // IMC transport step with Random Walk
    void transport(const MT &mesh, 
		   const Opacity<MT,Gray_Frequency> &opacity, 
		   Tally<MT> &tally, 
		   SP_Random_Walk = SP_Random_Walk(),
		   SP_Surface_tracker = SP_Surface_tracker(),
		   SP_Diagnostic = SP_Diagnostic()); 

    // >>> DIAGNOSTIC FUNCTIONS

    // Overloaded equals operator.
    bool operator==(const Gray_Particle<MT> &) const;

    //! Overloaded not equals operator.
    bool operator!=(const Gray_Particle<MT> &) const;
 
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
    Check (Base::r.size() < 4);
    Check (Base::omega.size() == 3);
    Check (Base::cell > 0); 
}

} // end namespace rtt_imc

#endif                          // RTT_imc_Gray_Particle_HH

//---------------------------------------------------------------------------//
//                              end of imc/Gray_Particle.hh
//---------------------------------------------------------------------------//
