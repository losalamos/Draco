//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Hybrid_Diffusion.hh
 * \author Thomas M. Evans
 * \date   Fri Feb 21 11:17:04 2003
 * \brief  Hybrid_Diffusion class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_imc_Hybrid_Diffusion_HH
#define RTT_imc_Hybrid_Diffusion_HH

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Hybrid_Diffusion
 *
 * \brief Base class for hybrid diffusion/IMC methods.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Hybrid_Diffusion 
{
  public:
    //! Hybrid diffusion/IMC methods.
    enum Methods {TRANSPORT   = 0,
		  RANDOM_WALK = 1,
		  DDIMC       = 2};

  public:
    //! Virtual Destructor.
    virtual inline ~Hybrid_Diffusion() = 0;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*! 
 * The virtual destructor must be defined, even though it is pure virtual.
 */
Hybrid_Diffusion::~Hybrid_Diffusion()
{
}

} // end namespace rtt_imc

#endif                          // RTT_imc_Hybrid_Diffusion_HH

//---------------------------------------------------------------------------//
//                              end of imc/Hybrid_Diffusion.hh
//---------------------------------------------------------------------------//
