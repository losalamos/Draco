//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   P1Diffusion/P1Momentum.hh
 * \author Randy M. Roberts
 * \date   Wed Feb  2 09:26:00 2000
 * \brief  Handles the velocity dependent parts of P1Diffusion
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __P1Diffusion_P1Momentum_hh__
#define __P1Diffusion_P1Momentum_hh__

namespace rtt_P1Diffusion
{
 
//===========================================================================//
/*!
 * \class P1Momentum
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, bool HASVELOCITY=true>
class P1Momentum 
{

    // NESTED CLASSES AND TYPEDEFS

  private:
     
    typedef typename MT::fcdsf fcdsf;
    typedef typename MT::ccsf ccsf;
    typedef typename MT::bssf bssf;
    typedef typename MT::ncvsf ncvsf;
    typedef typename MT::vcvsf vcvsf;
    typedef typename MT::vcsf vcsf;
    typedef typename MT::FieldConstructor FieldConstructor;

  public:

    typedef MT MeshType;
    typedef ccsf  IntensityField;
    typedef fcdsf FluxField;
    typedef fcdsf DiscFluxField;
    typedef ncvsf MomentumField;
    typedef vcvsf DiscMomentumField;
    typedef vcsf DiscKineticEnergyField;
    typedef fcdsf DiffCoefField;

    // DATA
    
    FieldConstructor fCtor;

  public:

    // CREATORS
    
    P1Momentum(const FieldConstructor &fCtor_)
	: fCtor(fCtor_)
    {
	// empty
    }
    
    //Defaulted: P1Momentum(const P1Momentum &rhs);
    //Defaulted: ~P1Momentum();

    // MANIPULATORS
    
    //Defaulted: P1Momentum& operator=(const P1Momentum &rhs);

    // ACCESSORS

    void discFluxToDiscMomentum(DiscMomentumField &result,
				const DiscFluxField &flux) const;

    void dotProduct(DiscKineticEnergyField &result,
		    const DiscMomentumField &vec1,
		    const DiscMomentumField &vec2) const;
     
    void dotProduct(DiscKineticEnergyField &KEnergy,
		    const DiscFluxField &sigmaF,
		    const DiscMomentumField &velocity) const;

  private:
    
    // IMPLEMENTATION
};

template<class MT>
class P1Momentum<MT, false>
{

    // NESTED CLASSES AND TYPEDEFS

  private:
     
    struct Dummy { /* empty */ };
    
    typedef typename MT::fcdsf fcdsf;
    typedef typename MT::ccsf ccsf;
    typedef typename MT::bssf bssf;
    typedef Dummy ncvsf;
    typedef Dummy vcvsf;
    typedef Dummy vcsf;
    typedef typename MT::FieldConstructor FieldConstructor;

  public:

    typedef MT MeshType;
    typedef ccsf  IntensityField;
    typedef fcdsf FluxField;
    typedef fcdsf DiscFluxField;
    typedef ncvsf MomentumField;
    typedef vcvsf DiscMomentumField;
    typedef vcsf DiscKineticEnergyField;
    typedef fcdsf DiffCoefField;

    // DATA
    
    FieldConstructor fCtor;

  public:

    // CREATORS
    
    P1Momentum(const FieldConstructor &fCtor_)
	: fCtor(fCtor_)
    {
	// empty
    }
    
    //Defaulted: P1Momentum(const P1Momentum &rhs);
    //Defaulted: ~P1Momentum();

    // MANIPULATORS
    
    //Defaulted: P1Momentum& operator=(const P1Momentum &rhs);

    // ACCESSORS

    void discFluxToDiscMomentum(DiscMomentumField &result,
				const DiscFluxField &flux) const
    {
	// empty
    }

    void dotProduct(DiscKineticEnergyField &result,
		    const DiscMomentumField &vec1,
		    const DiscMomentumField &vec2) const
    {
	// empty
    }
     
    void dotProduct(DiscKineticEnergyField &KEnergy,
		    const DiscFluxField &sigmaF,
		    const DiscMomentumField &velocity) const
    {
	// empty
    }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_P1Diffusion

#endif                          // __P1Diffusion_P1Momentum_hh__

//---------------------------------------------------------------------------//
//                              end of P1Diffusion/P1Momentum.hh
//---------------------------------------------------------------------------//
