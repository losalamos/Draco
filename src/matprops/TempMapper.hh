//----------------------------------*-C++-*----------------------------------//
// TempMapper.hh
// Randy M. Roberts
// Thu Oct  8 09:44:09 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_TempMapper_hh__
#define __matprops_TempMapper_hh__

#include "ds++/SP.hh"

namespace rtt_matprops
{

 // using is inside restricted namespace
 
 using dsxx::SP;
 
 //===========================================================================//
 // class TempMapper - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 template <class MT>
 class TempMapper
 {

     // NESTED CLASSES AND TYPEDEFS

     typedef typename MT::FieldConstructor FieldConstructor;

     // DATA

     double gamma;
     const SP<MT> &spMesh;
     const FieldConstructor &FC;
    
   public:

     // CREATORS
    
     TempMapper(const SP<MT> &spMesh_, double gamma_, const FieldConstructor &FC_);

     // MANIPULATORS
    
     // ACCESSORS
     
     template<class FCV, class CCV, class OpAssignV>
     void tempCC2FC(FCV &faceTempsByMat, const CCV &cellTempsByMat,
		    const typename MT::ccsf &cellTempsAvg,
		    const OpAssignV &opAssignV) const;

   private:
    
     // IMPLEMENTATION

     void calculateAlpha(typename MT::fcdsf &alpha,
			 const typename MT::fcdsf &faceTempsAvg) const;
 };

} // end namespace rtt_matprops

#endif                          // __matprops_TempMapper_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/TempMapper.hh
//---------------------------------------------------------------------------//
