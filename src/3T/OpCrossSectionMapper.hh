//----------------------------------*-C++-*----------------------------------//
// OpCrossSectionMapper.hh
// Randy M. Roberts
// Fri Nov 13 14:36:46 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_OpCrossSectionMapper_hh__
#define __3T_OpCrossSectionMapper_hh__

#include "3T/CrossSectionMapper.hh"

// DEFINING NAMESPACE

namespace rtt_3T {

 //===========================================================================//
 // class OpCrossSectionMapper - 

 // 
 //===========================================================================//

 template<class DS, class Op>
 class OpCrossSectionMapper : public CrossSectionMapper<DS>
 {

     // NESTED CLASSES AND TYPEDEFS

   private:

     typedef typename DS::MeshType MT;
#ifdef P13T_MOMENTUM_DEPOSITION
     typedef typename
     CrossSectionMapper<DS>::DiscKineticEnergyField DiscKineticEnergyField;
#endif
     typedef typename CrossSectionMapper<DS>::fcdsf fcdsf;

     // CREATORS

   public:
     
     // MANIPULATORS

     // ACCESSORS

#ifdef P13T_MOMENTUM_DEPOSITION
     virtual void mapCrossSections(DiscKineticEnergyField &vcSigma,
				   const fcdsf &fcSigma) const
     {
	 MT::scatter(vcSigma, fcSigma, Op());
     }
#endif
 };

} // end namespace rtt_3T

#endif                          // __3T_OpCrossSectionMapper_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/OpCrossSectionMapper.hh
//---------------------------------------------------------------------------//
