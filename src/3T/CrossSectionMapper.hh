//----------------------------------*-C++-*----------------------------------//
// CrossSectionMapper.hh
// Randy M. Roberts
// Fri Nov 13 14:36:46 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_CrossSectionMapper_hh__
#define __3T_CrossSectionMapper_hh__

#ifndef P13T_MOMENTUM_DEPOSITION
#define P13T_MOMENTUM_DEPOSITION
#endif

// DEFINING NAMESPACE

namespace rtt_3T {

 //===========================================================================//
 // class CrossSectionMapper - 

 // 
 //===========================================================================//

 template<class DS>
 class CrossSectionMapper
 {

     // NESTED CLASSES AND TYPEDEFS

   public:
     
#ifdef P13T_MOMENTUM_DEPOSITION
     typedef typename DS::DiscKineticEnergyField DiscKineticEnergyField;
#endif
     typedef typename DS::MeshType::fcdsf fcdsf;

     // CREATORS

   public:
     
     virtual ~CrossSectionMapper()
     {
	 // empty
     }
     
     // MANIPULATORS

     // ACCESSORS

#ifdef P13T_MOMENTUM_DEPOSITION
     virtual void mapCrossSections(DiscKineticEnergyField &vcSigma,
				   const fcdsf &fcSigma) const = 0;
#endif
 };

} // end namespace rtt_3T

#endif                          // __3T_CrossSectionMapper_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/CrossSectionMapper.hh
//---------------------------------------------------------------------------//
