//----------------------------------*-C++-*----------------------------------//
// CrossSectionMapper.hh
// Randy M. Roberts
// Fri Nov 13 14:36:46 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_CrossSectionMapper_hh__
#define __3T_CrossSectionMapper_hh__

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
     
     typedef typename DS::DiscKineticEnergyField DiscKineticEnergyField;
     typedef typename DS::MeshType::fcdsf fcdsf;

     // CREATORS

   public:
     
     virtual ~CrossSectionMapper()
     {
	 // empty
     }
     
     // MANIPULATORS

     // ACCESSORS

     virtual void mapCrossSections(DiscKineticEnergyField &vcSigma,
				   const fcdsf &fcSigma) const = 0;
 };

} // end namespace rtt_3T

#endif                          // __3T_CrossSectionMapper_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/CrossSectionMapper.hh
//---------------------------------------------------------------------------//
